#include "dfs_solver.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>


static std::uint64_t splitmix64(std::uint64_t x) {
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

static std::string to_lower_ascii(const std::string& input) {
	std::string out;
	out.reserve(input.size());
	for (char ch : input) {
		out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
	}
	return out;
}

const char* to_string(memo_mode mode) {
	switch (mode) {
	case memo_mode::none:
		return "none";
	case memo_mode::passive:
		return "passive";
	case memo_mode::predictive:
		return "predictive";
	case memo_mode::solution:
		return "solution";
	default:
		return "unknown";
	}
}

const char* to_string(memo_cleaning_policy policy) {
	switch (policy) {
	case memo_cleaning_policy::lru:
		return "lru";
	case memo_cleaning_policy::lufo:
		return "lufo";
	default:
		return "unknown";
	}
}

const char* to_string(branch_rule rule) {
	switch (rule) {
	case branch_rule::slack:
		return "slack";
	case branch_rule::edd:
		return "edd";
	case branch_rule::spt:
		return "spt";
	case branch_rule::decomposition:
		return "decomposition";
	default:
		return "unknown";
	}
}

bool parse_memo_mode(const std::string& text, memo_mode& out) {
	const std::string key = to_lower_ascii(text);
	if (key == "none") {
		out = memo_mode::none;
		return true;
	}
	if (key == "passive") {
		out = memo_mode::passive;
		return true;
	}
	if (key == "predictive") {
		out = memo_mode::predictive;
		return true;
	}
	if (key == "solution") {
		out = memo_mode::solution;
		return true;
	}
	return false;
}

bool parse_memo_cleaning_policy(const std::string& text, memo_cleaning_policy& out) {
	const std::string key = to_lower_ascii(text);
	if (key == "lru") {
		out = memo_cleaning_policy::lru;
		return true;
	}
	if (key == "lufo") {
		out = memo_cleaning_policy::lufo;
		return true;
	}
	return false;
}

bool parse_branch_rule(const std::string& text, branch_rule& out) {
	const std::string key = to_lower_ascii(text);
	if (key == "slack") {
		out = branch_rule::slack;
		return true;
	}
	if (key == "edd") {
		out = branch_rule::edd;
		return true;
	}
	if (key == "spt") {
		out = branch_rule::spt;
		return true;
	}
	if (key == "decomposition") {
		out = branch_rule::decomposition;
		return true;
	}
	return false;
}

dfs_solver::dfs_solver(dfs_config config)
	: config_(config), memo_(config.memo_capacity, config.memo_cleaning) {}

solve_result dfs_solver::solve(const instance& inst) {
	std::string error;
	if (!validate_instance(inst, &error)) {
		throw std::invalid_argument(error);
	}

	inst_ = &inst;
	n_ = static_cast<int>(inst.jobs.size());

	perm_jobs_.resize(n_);
	std::iota(perm_jobs_.begin(), perm_jobs_.end(), 0);

	const int words = (n_ + 63) / 64;
	remaining_bits_.assign(static_cast<std::size_t>(words), 0ULL);
	for (int j = 0; j < n_; ++j) {
		remaining_bits_[static_cast<std::size_t>(j >> 6)] |= (1ULL << (j & 63));
	}

	zobrist_job_.resize(n_);
	subset_hash_ = 0ULL;
	std::uint64_t seed = config_.zobrist_seed;
	for (int j = 0; j < n_; ++j) {
		seed = splitmix64(seed + static_cast<std::uint64_t>(j + 1));
		zobrist_job_[j] = seed;
		subset_hash_ ^= seed;
	}

	memo_.set_cleaning_policy(config_.memo_cleaning, false);
	memo_.clear();
	memo_.set_capacity(config_.memo_capacity, false);
	stats_ = {};

	const bool use_exact_memo_main =
		(config_.memo == memo_mode::predictive) || (config_.memo == memo_mode::solution);
	const bool use_lb_memo_main = (config_.memo == memo_mode::predictive);
	const auto t0 = std::chrono::steady_clock::now();
	const long long best_cost = solve_state(0, 0, nullptr, use_exact_memo_main, use_lb_memo_main, true);
	const auto t1 = std::chrono::steady_clock::now();

	stats_.elapsed_ms =
		std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t1 - t0).count();

	const memo_table_stats memo_stats = memo_.stats();
	stats_.memo_hits = memo_stats.hits;
	stats_.memo_misses = memo_stats.misses;
	stats_.memo_inserts = memo_stats.inserts;
	stats_.memo_updates = memo_stats.updates;
	stats_.memo_evictions = memo_stats.evictions;
	stats_.memo_peak_size = memo_stats.peak_size;
	stats_.memo_final_size = memo_stats.final_size;

	schedule best;
	best.cost = best_cost;
	if (config_.reconstruct_order && n_ <= 200) {
		best.order = reconstruct_order(best_cost);
	}
	if (!best.order.empty() && evaluate_sum_tardiness(inst, best.order) != best.cost) {
		best.order.clear();
	}

	return {best, stats_};
}

long long dfs_solver::solve_state(int depth, int current_time, const memo_lookup_result* known_lookup,
	bool use_exact_memo, bool use_lb_memo, bool track_stats) {
	if (track_stats) {
		++stats_.nodes;
	}

	if (depth == n_) {
		if (track_stats) {
			++stats_.leaves;
		}
		return 0;
	}

	const bool memo_available = (config_.memo != memo_mode::none) && (config_.memo_capacity > 0);
	memo_lookup_result state_lookup;
	bool has_lookup = false;

	if (use_exact_memo && memo_available) {
		if (known_lookup != nullptr) {
			state_lookup = *known_lookup;
			has_lookup = state_lookup.found;
		}
		else {
			state_lookup = memo_.lookup(remaining_bits_, current_time, subset_hash_, track_stats);
			has_lookup = state_lookup.found;
		}
		if (has_lookup && state_lookup.has_exact) {
			if (track_stats) {
				++stats_.pruned_by_memo_exact;
			}
			return state_lookup.exact;
		}
	}

	const long long base_lb = lower_bound_additional(depth, current_time);
	long long effective_lb = base_lb;
	if (use_lb_memo && has_lookup && state_lookup.lower_bound > effective_lb) {
		effective_lb = state_lookup.lower_bound;
	}

	if (use_lb_memo && memo_available) {
		memo_.store_lower_bound(remaining_bits_, current_time, subset_hash_, effective_lb, track_stats);
	}

	int heuristic_first_job = -1;
	long long best_cost = heuristic_upper_bound_edd(depth, current_time, &heuristic_first_job);
	int best_first_job = heuristic_first_job;

	if (effective_lb >= best_cost) {
		if (memo_available && use_exact_memo) {
			memo_.store_exact(remaining_bits_, current_time, subset_hash_, best_cost, best_first_job, track_stats);
		}
		return best_cost;
	}

	const std::vector<int> candidates = build_branch_indices(depth, current_time);
	for (int idx : candidates) {
		const int job_idx = perm_jobs_[idx];
		std::swap(perm_jobs_[depth], perm_jobs_[idx]);

		const job& j = inst_->jobs[job_idx];
		const int child_time = current_time + j.p;
		const long long immediate = tardiness(child_time, j.d);

		remove_job_from_state(job_idx);

		const long long child_lb_base = lower_bound_additional(depth + 1, child_time);
		long long child_lb = immediate + child_lb_base;

		memo_lookup_result child_lookup;
		bool child_lookup_found = false;
		bool memo_lb_tightened = false;

		if (use_exact_memo && memo_available) {
			child_lookup = memo_.lookup(remaining_bits_, child_time, subset_hash_, track_stats);
			child_lookup_found = child_lookup.found;

			if (child_lookup_found && child_lookup.has_exact) {
				const long long candidate = immediate + child_lookup.exact;
				if (candidate < best_cost) {
					best_cost = candidate;
					best_first_job = job_idx;
				}
				if (track_stats) {
					++stats_.pruned_by_memo_exact;
				}
				restore_job_to_state(job_idx);
				std::swap(perm_jobs_[depth], perm_jobs_[idx]);
				continue;
			}

			if (use_lb_memo && child_lookup_found && child_lookup.lower_bound > child_lb_base) {
				child_lb = immediate + child_lookup.lower_bound;
				memo_lb_tightened = true;
			}
		}

		if (child_lb >= best_cost) {
			if (track_stats) {
				++stats_.pruned_by_bound;
				if (memo_lb_tightened) {
					++stats_.pruned_by_memo_lb;
				}
			}
			restore_job_to_state(job_idx);
			std::swap(perm_jobs_[depth], perm_jobs_[idx]);
			continue;
		}

		const long long child_exact = solve_state(depth + 1, child_time,
			child_lookup_found ? &child_lookup : nullptr, use_exact_memo, use_lb_memo, track_stats);
		const long long candidate = immediate + child_exact;
		if (candidate < best_cost) {
			best_cost = candidate;
			best_first_job = job_idx;
		}

		restore_job_to_state(job_idx);
		std::swap(perm_jobs_[depth], perm_jobs_[idx]);
	}

	if (memo_available && use_exact_memo) {
		memo_.store_exact(remaining_bits_, current_time, subset_hash_, best_cost, best_first_job, track_stats);
	}

	return best_cost;
}

long long dfs_solver::lower_bound_additional(int depth, int current_time) const {
	long long lb = 0;
	for (int i = depth; i < n_; ++i) {
		const job& j = inst_->jobs[perm_jobs_[i]];
		lb += tardiness(current_time + j.p, j.d);
	}
	return lb;
}

long long dfs_solver::heuristic_upper_bound_edd(int depth, int current_time, int* first_job) const {
	std::vector<int> jobs;
	jobs.reserve(static_cast<std::size_t>(n_ - depth));
	for (int i = depth; i < n_; ++i) {
		jobs.push_back(perm_jobs_[i]);
	}

	std::sort(jobs.begin(), jobs.end(), [this](int a, int b) {
		const job& ja = inst_->jobs[a];
		const job& jb = inst_->jobs[b];
		if (ja.d != jb.d) {
			return ja.d < jb.d;
		}
		if (ja.p != jb.p) {
			return ja.p < jb.p;
		}
		return a < b;
	});

	if (first_job != nullptr) {
		*first_job = jobs.empty() ? -1 : jobs.front();
	}

	int t = current_time;
	long long sum = 0;
	for (int j : jobs) {
		t += inst_->jobs[j].p;
		sum += tardiness(t, inst_->jobs[j].d);
	}
	return sum;
}

std::vector<int> dfs_solver::build_branch_indices(int depth, int current_time) const {
	std::vector<int> idx;
	idx.reserve(static_cast<std::size_t>(n_ - depth));
	for (int i = depth; i < n_; ++i) {
		idx.push_back(i);
	}

	switch (config_.branching) {
	case branch_rule::slack:
		std::sort(idx.begin(), idx.end(), [this, current_time](int ia, int ib) {
			const job& a = inst_->jobs[perm_jobs_[ia]];
			const job& b = inst_->jobs[perm_jobs_[ib]];
			const long long slack_a = static_cast<long long>(a.d) - (current_time + a.p);
			const long long slack_b = static_cast<long long>(b.d) - (current_time + b.p);
			if (slack_a != slack_b) {
				return slack_a < slack_b;
			}
			if (a.d != b.d) {
				return a.d < b.d;
			}
			if (a.p != b.p) {
				return a.p > b.p;
			}
			return perm_jobs_[ia] < perm_jobs_[ib];
		});
		break;
	case branch_rule::edd:
		std::sort(idx.begin(), idx.end(), [this](int ia, int ib) {
			const job& a = inst_->jobs[perm_jobs_[ia]];
			const job& b = inst_->jobs[perm_jobs_[ib]];
			if (a.d != b.d) {
				return a.d < b.d;
			}
			if (a.p != b.p) {
				return a.p < b.p;
			}
			return perm_jobs_[ia] < perm_jobs_[ib];
		});
		break;
	case branch_rule::spt:
		std::sort(idx.begin(), idx.end(), [this](int ia, int ib) {
			const job& a = inst_->jobs[perm_jobs_[ia]];
			const job& b = inst_->jobs[perm_jobs_[ib]];
			if (a.p != b.p) {
				return a.p < b.p;
			}
			if (a.d != b.d) {
				return a.d < b.d;
			}
			return perm_jobs_[ia] < perm_jobs_[ib];
		});
		break;
	case branch_rule::decomposition:
		std::sort(idx.begin(), idx.end(), [this, current_time](int ia, int ib) {
			const job& a = inst_->jobs[perm_jobs_[ia]];
			const job& b = inst_->jobs[perm_jobs_[ib]];
			const long long key_a = static_cast<long long>(a.d) - (current_time + a.p);
			const long long key_b = static_cast<long long>(b.d) - (current_time + b.p);
			if (key_a != key_b) {
				return key_a < key_b;
			}
			if (a.p != b.p) {
				return a.p > b.p;
			}
			if (a.d != b.d) {
				return a.d < b.d;
			}
			return perm_jobs_[ia] < perm_jobs_[ib];
		});
		break;
	}

	return idx;
}

std::vector<int> dfs_solver::reconstruct_order(long long optimal_cost) {
	std::vector<int> order;
	order.reserve(static_cast<std::size_t>(n_));

	const std::vector<int> saved_perm = perm_jobs_;
	const std::vector<std::uint64_t> saved_bits = remaining_bits_;
	const std::uint64_t saved_hash = subset_hash_;

	int current_time = 0;
	for (int depth = 0; depth < n_; ++depth) {
		int chosen_job = -1;

		if (config_.memo_capacity > 0) {
			const memo_lookup_result lookup = memo_.lookup(remaining_bits_, current_time, subset_hash_, false);
			if (lookup.found && lookup.has_exact && lookup.best_job >= 0 && is_job_remaining(lookup.best_job)) {
				chosen_job = lookup.best_job;
			}
		}

		if (chosen_job < 0) {
			long long best = std::numeric_limits<long long>::max() / 4;
			const std::vector<int> candidates = build_branch_indices(depth, current_time);
			for (int idx : candidates) {
				const int jidx = perm_jobs_[idx];
				std::swap(perm_jobs_[depth], perm_jobs_[idx]);

				const job& j = inst_->jobs[jidx];
				const int child_time = current_time + j.p;
				const long long immediate = tardiness(child_time, j.d);

				remove_job_from_state(jidx);

				memo_lookup_result child_lookup;
				bool child_lookup_found = false;
				if (config_.memo_capacity > 0) {
					child_lookup = memo_.lookup(remaining_bits_, child_time, subset_hash_, false);
					child_lookup_found = child_lookup.found;
				}

				long long child_exact = 0;
				if (child_lookup_found && child_lookup.has_exact) {
					child_exact = child_lookup.exact;
				}
				else {
					const bool use_exact_memo = (config_.memo == memo_mode::predictive) ||
						(config_.memo == memo_mode::solution);
					const bool use_lb_memo = (config_.memo == memo_mode::predictive);
					child_exact = solve_state(depth + 1, child_time, child_lookup_found ? &child_lookup : nullptr,
						use_exact_memo, use_lb_memo, false);
				}

				const long long candidate = immediate + child_exact;
				if (candidate < best) {
					best = candidate;
					chosen_job = jidx;
				}

				restore_job_to_state(jidx);
				std::swap(perm_jobs_[depth], perm_jobs_[idx]);
			}
		}

		if (chosen_job < 0) {
			order.clear();
			break;
		}

		int chosen_pos = -1;
		for (int pos = depth; pos < n_; ++pos) {
			if (perm_jobs_[pos] == chosen_job) {
				chosen_pos = pos;
				break;
			}
		}
		if (chosen_pos < 0) {
			order.clear();
			break;
		}

		std::swap(perm_jobs_[depth], perm_jobs_[chosen_pos]);
		remove_job_from_state(chosen_job);
		order.push_back(chosen_job);
		current_time += inst_->jobs[chosen_job].p;
	}

	perm_jobs_ = saved_perm;
	remaining_bits_ = saved_bits;
	subset_hash_ = saved_hash;

	if (order.size() != static_cast<std::size_t>(n_)) {
		order.clear();
	}
	if (!order.empty() && evaluate_sum_tardiness(*inst_, order) != static_cast<schedule_cost_t>(optimal_cost)) {
		order.clear();
	}

	return order;
}

void dfs_solver::remove_job_from_state(int job_idx) {
	const std::size_t word = static_cast<std::size_t>(job_idx >> 6);
	const std::uint64_t mask = (1ULL << (job_idx & 63));
	remaining_bits_[word] &= ~mask;
	subset_hash_ ^= zobrist_job_[job_idx];
}

void dfs_solver::restore_job_to_state(int job_idx) {
	const std::size_t word = static_cast<std::size_t>(job_idx >> 6);
	const std::uint64_t mask = (1ULL << (job_idx & 63));
	remaining_bits_[word] |= mask;
	subset_hash_ ^= zobrist_job_[job_idx];
}

bool dfs_solver::is_job_remaining(int job_idx) const {
	const std::size_t word = static_cast<std::size_t>(job_idx >> 6);
	const std::uint64_t mask = (1ULL << (job_idx & 63));
	return (remaining_bits_[word] & mask) != 0ULL;
}
