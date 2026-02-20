#include <algorithm>
#include <cctype>
#include <cmath>
#include <csignal>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "core.h"
#include "dfs_solver.h"
#include "generator.h"

namespace {
volatile std::sig_atomic_t g_stop_requested = 0;

void on_interrupt_signal(int) {
	g_stop_requested = 1;
}

enum class due_date_model {
	uniform = 0,
	potts = 1
};

const char* to_string(due_date_model model) {
	switch (model) {
	case due_date_model::uniform:
		return "uniform";
	case due_date_model::potts:
		return "potts";
	default:
		return "unknown";
	}
}

bool parse_due_date_model(const std::string& text, due_date_model& out) {
	std::string key = text;
	for (char& ch : key) {
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	if (key == "uniform") {
		out = due_date_model::uniform;
		return true;
	}
	if (key == "potts") {
		out = due_date_model::potts;
		return true;
	}
	return false;
}

struct cli_options {
	bool show_help = false;
	bool self_check = false;
	bool benchmark = false;
	bool bench_all_configs = false;
	bool memo_cap_explicit = false;
	std::string file_path;
	std::string bench_csv_path;
	int n = 16;
	int n_from = 16;
	int n_to = 16;
	int n_step = 1;
	int instances = 10;
	bool n_from_explicit = false;
	bool n_to_explicit = false;
	int p_min = 1;
	int p_max = 30;
	int d_min = 1;
	int d_max = -1;
	bool d_max_explicit = false;
	due_date_model due_model = due_date_model::potts;
	double due_range = 0.6;
	double due_tardiness = 0.6;
	std::size_t memory_budget_mb = 1024;
	dfs_config config{};
};

std::size_t state_count_upper_bound(std::size_t n) {
	if (n >= std::numeric_limits<std::size_t>::digits) {
		return std::numeric_limits<std::size_t>::max();
	}
	return (std::size_t{1} << n);
}

std::size_t estimate_memo_entry_bytes(std::size_t n) {
	const std::size_t words = (n + 63) / 64;
	const std::size_t bits_payload = words * sizeof(std::uint64_t);
	// Approximate per-entry footprint with container/node and index overhead.
	// The previous estimate was too optimistic and could underutilize or overshoot RAM.
	return 320 + bits_payload;
}

std::size_t memo_capacity_for_budget(std::size_t n, std::size_t memory_budget_mb) {
	if (memory_budget_mb == 0) {
		return 0;
	}

	const std::size_t mb_to_bytes = std::size_t{1024} * std::size_t{1024};
	const std::size_t max_mb = std::numeric_limits<std::size_t>::max() / mb_to_bytes;
	const std::size_t safe_mb = std::min(memory_budget_mb, max_mb);
	const std::size_t budget_bytes = safe_mb * mb_to_bytes;
	const std::size_t bytes_per_entry = std::max<std::size_t>(estimate_memo_entry_bytes(n), 1);
	const std::size_t usable_bytes = static_cast<std::size_t>(budget_bytes * 0.60);
	std::size_t capacity = usable_bytes / bytes_per_entry;
	if (capacity == 0) {
		capacity = 1;
	}

	const std::size_t full_state_count = state_count_upper_bound(n);
	if (capacity > full_state_count) {
		capacity = full_state_count;
	}
	return capacity;
}

std::size_t recommended_memo_capacity(std::size_t n, std::size_t memory_budget_mb) {
	if (n == 0) {
		return 1;
	}
	return memo_capacity_for_budget(n, memory_budget_mb);
}

int default_due_date_upper(int n, int p_max, int d_min) {
	const long long scaled = (static_cast<long long>(n) * static_cast<long long>(p_max)) / 2LL;
	const long long capped = std::min<long long>(scaled, static_cast<long long>(std::numeric_limits<int>::max()));
	return static_cast<int>(std::max<long long>(d_min, capped));
}

void print_usage() {
	std::cout << "Usage:\n"
		<< "  kursovaya.exe [--file <path>] [--memo none|passive|predictive|solution]\n"
		<< "                [--memo-cap <int>] [--memo-clean lru|lufo] [--branch slack|edd|spt|decomposition] [--seed <uint64>]\n"
		<< "                [--n <int>] [--p-min <int>] [--p-max <int>] [--d-min <int>] [--d-max <int>]\n"
		<< "                [--due-model uniform|potts] [--due-range <double>] [--due-tardiness <double>]\n"
		<< "                [--mem-budget-mb <int>] [--no-reconstruct] [--self-check]\n"
		<< "                [--bench-csv <path>] [--n-from <int>] [--n-to <int>] [--n-step <int>]\n"
		<< "                [--instances <int>] [--bench-all-configs] [--help]\n\n"
		<< "Input file format:\n"
		<< "  n\n"
		<< "  p_0 d_0\n"
		<< "  ...\n"
		<< "  p_{n-1} d_{n-1}\n\n"
		<< "Benchmark mode:\n"
		<< "  --bench-csv enables batch runs and appends one CSV row per solved instance.\n"
		<< "  Rows are flushed immediately, so partial results survive interruption.\n";
}

bool parse_int_arg(const char* text, int& out) {
	try {
		const int v = std::stoi(text);
		out = v;
		return true;
	}
	catch (...) {
		return false;
	}
}

bool parse_u64_arg(const char* text, std::uint64_t& out) {
	try {
		const std::uint64_t v = static_cast<std::uint64_t>(std::stoull(text));
		out = v;
		return true;
	}
	catch (...) {
		return false;
	}
}

bool parse_size_arg(const char* text, std::size_t& out) {
	try {
		const std::size_t v = static_cast<std::size_t>(std::stoull(text));
		out = v;
		return true;
	}
	catch (...) {
		return false;
	}
}

bool parse_double_arg(const char* text, double& out) {
	try {
		const double v = std::stod(text);
		out = v;
		return true;
	}
	catch (...) {
		return false;
	}
}

bool parse_cli(int argc, char** argv, cli_options& opts, std::string& error) {
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		auto require_value = [&](const char* name) -> const char* {
			if (i + 1 >= argc) {
				error = std::string("Missing value for ") + name + ".";
				return nullptr;
			}
			++i;
			return argv[i];
		};

		if (arg == "--help" || arg == "-h") {
			opts.show_help = true;
		}
		else if (arg == "--self-check") {
			opts.self_check = true;
		}
		else if (arg == "--no-reconstruct") {
			opts.config.reconstruct_order = false;
		}
		else if (arg == "--file") {
			const char* value = require_value("--file");
			if (value == nullptr) {
				return false;
			}
			opts.file_path = value;
		}
		else if (arg == "--memo") {
			const char* value = require_value("--memo");
			if (value == nullptr) {
				return false;
			}
			memo_mode mode{};
			if (!parse_memo_mode(value, mode)) {
				error = "Unknown memo mode: " + std::string(value);
				return false;
			}
			opts.config.memo = mode;
		}
		else if (arg == "--branch") {
			const char* value = require_value("--branch");
			if (value == nullptr) {
				return false;
			}
			branch_rule rule{};
			if (!parse_branch_rule(value, rule)) {
				error = "Unknown branch rule: " + std::string(value);
				return false;
			}
			opts.config.branching = rule;
		}
		else if (arg == "--memo-cap") {
			const char* value = require_value("--memo-cap");
			if (value == nullptr) {
				return false;
			}
			if (!parse_size_arg(value, opts.config.memo_capacity)) {
				error = "Invalid --memo-cap value.";
				return false;
			}
			opts.memo_cap_explicit = true;
		}
		else if (arg == "--memo-clean") {
			const char* value = require_value("--memo-clean");
			if (value == nullptr) {
				return false;
			}
			memo_cleaning_policy policy{};
			if (!parse_memo_cleaning_policy(value, policy)) {
				error = "Unknown --memo-clean value: " + std::string(value);
				return false;
			}
			opts.config.memo_cleaning = policy;
		}
		else if (arg == "--mem-budget-mb") {
			const char* value = require_value("--mem-budget-mb");
			if (value == nullptr) {
				return false;
			}
			if (!parse_size_arg(value, opts.memory_budget_mb)) {
				error = "Invalid --mem-budget-mb value.";
				return false;
			}
		}
		else if (arg == "--seed") {
			const char* value = require_value("--seed");
			if (value == nullptr) {
				return false;
			}
			if (!parse_u64_arg(value, opts.config.zobrist_seed)) {
				error = "Invalid --seed value.";
				return false;
			}
		}
		else if (arg == "--n") {
			const char* value = require_value("--n");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.n)) {
				error = "Invalid --n value.";
				return false;
			}
		}
		else if (arg == "--n-from") {
			const char* value = require_value("--n-from");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.n_from)) {
				error = "Invalid --n-from value.";
				return false;
			}
			opts.n_from_explicit = true;
		}
		else if (arg == "--n-to") {
			const char* value = require_value("--n-to");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.n_to)) {
				error = "Invalid --n-to value.";
				return false;
			}
			opts.n_to_explicit = true;
		}
		else if (arg == "--n-step") {
			const char* value = require_value("--n-step");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.n_step)) {
				error = "Invalid --n-step value.";
				return false;
			}
		}
		else if (arg == "--instances") {
			const char* value = require_value("--instances");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.instances)) {
				error = "Invalid --instances value.";
				return false;
			}
		}
		else if (arg == "--p-min") {
			const char* value = require_value("--p-min");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.p_min)) {
				error = "Invalid --p-min value.";
				return false;
			}
		}
		else if (arg == "--p-max") {
			const char* value = require_value("--p-max");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.p_max)) {
				error = "Invalid --p-max value.";
				return false;
			}
		}
		else if (arg == "--d-min") {
			const char* value = require_value("--d-min");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.d_min)) {
				error = "Invalid --d-min value.";
				return false;
			}
		}
		else if (arg == "--d-max") {
			const char* value = require_value("--d-max");
			if (value == nullptr) {
				return false;
			}
			if (!parse_int_arg(value, opts.d_max)) {
				error = "Invalid --d-max value.";
				return false;
			}
			opts.d_max_explicit = true;
		}
		else if (arg == "--due-model") {
			const char* value = require_value("--due-model");
			if (value == nullptr) {
				return false;
			}
			due_date_model model{};
			if (!parse_due_date_model(value, model)) {
				error = "Unknown --due-model value: " + std::string(value);
				return false;
			}
			opts.due_model = model;
		}
		else if (arg == "--due-range") {
			const char* value = require_value("--due-range");
			if (value == nullptr) {
				return false;
			}
			if (!parse_double_arg(value, opts.due_range)) {
				error = "Invalid --due-range value.";
				return false;
			}
		}
		else if (arg == "--due-tardiness") {
			const char* value = require_value("--due-tardiness");
			if (value == nullptr) {
				return false;
			}
			if (!parse_double_arg(value, opts.due_tardiness)) {
				error = "Invalid --due-tardiness value.";
				return false;
			}
		}
		else if (arg == "--bench-csv") {
			const char* value = require_value("--bench-csv");
			if (value == nullptr) {
				return false;
			}
			opts.benchmark = true;
			opts.bench_csv_path = value;
		}
		else if (arg == "--bench-all-configs") {
			opts.bench_all_configs = true;
		}
		else {
			error = "Unknown argument: " + arg;
			return false;
		}
	}

	if (opts.n < 0) {
		error = "--n must be non-negative.";
		return false;
	}
	if (opts.n_step <= 0) {
		error = "--n-step must be positive.";
		return false;
	}
	if (opts.instances <= 0) {
		error = "--instances must be positive.";
		return false;
	}

	if (opts.benchmark) {
		if (!opts.n_from_explicit && !opts.n_to_explicit) {
			opts.n_from = opts.n;
			opts.n_to = opts.n;
		}
		else if (!opts.n_from_explicit) {
			opts.n_from = opts.n_to;
		}
		else if (!opts.n_to_explicit) {
			opts.n_to = opts.n_from;
		}

		if (opts.n_from < 0 || opts.n_to < 0 || opts.n_to < opts.n_from) {
			error = "Invalid benchmark n-range.";
			return false;
		}
		if (opts.bench_csv_path.empty()) {
			error = "--bench-csv requires a file path.";
			return false;
		}
		if (!opts.file_path.empty()) {
			error = "--file cannot be used together with --bench-csv.";
			return false;
		}
	}
	else if (opts.bench_all_configs) {
		error = "--bench-all-configs requires --bench-csv.";
		return false;
	}

	if (opts.p_min <= 0 || opts.p_max < opts.p_min) {
		error = "Invalid processing-time range.";
		return false;
	}
	if (opts.p_max > static_cast<int>(std::numeric_limits<processing_time_t>::max())) {
		error = "Processing-time range exceeds supported type limits.";
		return false;
	}
	if (opts.due_model == due_date_model::uniform) {
		if (opts.d_min < 0) {
			error = "Due-date range cannot be negative.";
			return false;
		}
		if (opts.d_max < 0) {
			const int reference_n = opts.benchmark ? opts.n_to : opts.n;
			opts.d_max = default_due_date_upper(reference_n, opts.p_max, opts.d_min);
		}
		if (opts.d_max < opts.d_min) {
			error = "Invalid due-date range.";
			return false;
		}
		if (opts.d_max > static_cast<int>(std::numeric_limits<due_date_t>::max())) {
			error = "Due-date range exceeds supported type limits.";
			return false;
		}
	}
	else {
		if (!std::isfinite(opts.due_range) || opts.due_range <= 0.0) {
			error = "--due-range must be positive.";
			return false;
		}
		if (!std::isfinite(opts.due_tardiness) || opts.due_tardiness < 0.0 || opts.due_tardiness > 1.5) {
			error = "--due-tardiness must be in [0, 1.5].";
			return false;
		}
	}

	return true;
}

long long brute_force_optimal(const instance& inst, std::vector<int>* best_order = nullptr) {
	const int n = static_cast<int>(inst.jobs.size());
	std::vector<int> order(static_cast<std::size_t>(n));
	std::iota(order.begin(), order.end(), 0);

	long long best = std::numeric_limits<long long>::max();
	std::vector<int> best_local;

	do {
		const long long value = evaluate_sum_tardiness(inst, order);
		if (value < best) {
			best = value;
			best_local = order;
		}
	} while (std::next_permutation(order.begin(), order.end()));

	if (best_order != nullptr) {
		*best_order = std::move(best_local);
	}
	return best;
}

bool file_has_content(const std::string& path) {
	std::ifstream in(path, std::ios::binary | std::ios::ate);
	if (!in) {
		return false;
	}
	return in.tellg() > 0;
}

std::string csv_escape(std::string text) {
	bool needs_quotes = false;
	for (char ch : text) {
		if (ch == '"' || ch == ',' || ch == '\n' || ch == '\r') {
			needs_quotes = true;
			break;
		}
	}
	if (!needs_quotes) {
		return text;
	}

	std::string out;
	out.reserve(text.size() + 8);
	out.push_back('"');
	for (char ch : text) {
		if (ch == '"') {
			out.push_back('"');
		}
		out.push_back(ch);
	}
	out.push_back('"');
	return out;
}

instance generate_instance_for_run(const cli_options& opts, int n, std::uint64_t seed, int uniform_d_max_override) {
	if (opts.due_model == due_date_model::potts) {
		return generate_potts_instance(n, opts.p_min, opts.p_max, opts.due_range, opts.due_tardiness, seed);
	}

	int d_max_for_instance = uniform_d_max_override;
	if (!opts.d_max_explicit) {
		d_max_for_instance = default_due_date_upper(n, opts.p_max, opts.d_min);
	}
	return generate_random_instance(n, opts.p_min, opts.p_max, opts.d_min, d_max_for_instance, seed);
}

bool run_benchmark_csv(const cli_options& base, std::string& error) {
	std::ofstream csv(base.bench_csv_path, std::ios::out | std::ios::app);
	if (!csv) {
		error = "Failed to open CSV file: " + base.bench_csv_path;
		return false;
	}

	if (!file_has_content(base.bench_csv_path)) {
		csv << "n,instance_idx,seed,due_model,due_range,due_tardiness,memo_mode,memo_cleaning,branch_rule,memo_capacity,best_sum_tardiness,time_ms,"
			<< "nodes,leaves,pruned_by_bound,pruned_by_memo_exact,pruned_by_memo_lb,memo_hits,memo_misses,"
			<< "memo_inserts,memo_updates,memo_evictions,memo_peak_size,memo_final_size,status,error\n";
		csv.flush();
	}

	std::vector<memo_mode> modes;
	std::vector<branch_rule> branches;
	if (base.bench_all_configs) {
		modes = {memo_mode::none, memo_mode::passive, memo_mode::predictive, memo_mode::solution};
		branches = {branch_rule::slack, branch_rule::edd, branch_rule::spt, branch_rule::decomposition};
	}
	else {
		modes = {base.config.memo};
		branches = {base.config.branching};
	}

	g_stop_requested = 0;
	std::signal(SIGINT, on_interrupt_signal);

	std::cout << "[bench] Writing results to: " << base.bench_csv_path << "\n";
	std::cout << "[bench] Press Ctrl+C to stop safely after current case.\n";
#ifndef NDEBUG
	std::cout << "[bench] Warning: Debug build is active; benchmark timings can be much slower than Release.\n";
#endif

	for (int n = base.n_from; n <= base.n_to && !g_stop_requested; n += base.n_step) {
		for (int instance_idx = 0; instance_idx < base.instances && !g_stop_requested; ++instance_idx) {
			for (memo_mode mode : modes) {
				for (branch_rule branch : branches) {
					if (g_stop_requested) {
						break;
					}

					const std::uint64_t seed = base.config.zobrist_seed +
						static_cast<std::uint64_t>(n) * 1000003ULL +
						static_cast<std::uint64_t>(instance_idx) * 9176ULL +
						static_cast<std::uint64_t>(mode) * 271ULL +
						static_cast<std::uint64_t>(branch) * 37ULL;

					instance inst = generate_instance_for_run(base, n, seed, base.d_max);

					dfs_config cfg = base.config;
					cfg.memo = mode;
					cfg.branching = branch;
					cfg.reconstruct_order = false;

					if (cfg.memo == memo_mode::none) {
						cfg.memo_capacity = 0;
					}
					else if (!base.memo_cap_explicit) {
						cfg.memo_capacity = recommended_memo_capacity(inst.jobs.size(), base.memory_budget_mb);
					}

					schedule_cost_t cost = 0;
					solver_stats stats{};
					std::string status = "ok";
					std::string error_text;
					std::cout << "[bench] start n=" << n << " instance=" << instance_idx
						<< " due=" << to_string(base.due_model)
						<< " memo=" << to_string(mode) << " clean=" << to_string(cfg.memo_cleaning)
						<< " branch=" << to_string(branch) << "\n";
					std::cout.flush();

					try {
						dfs_solver solver(cfg);
						const solve_result result = solver.solve(inst);
						cost = result.best.cost;
						stats = result.stats;
					}
					catch (const std::exception& ex) {
						status = "error";
						error_text = ex.what();
					}

					csv << n << "," << instance_idx << "," << seed << ","
						<< to_string(base.due_model) << "," << base.due_range << "," << base.due_tardiness << ","
						<< to_string(mode) << ","
						<< to_string(cfg.memo_cleaning) << "," << to_string(branch) << "," << cfg.memo_capacity << "," << cost << ","
						<< stats.elapsed_ms << "," << stats.nodes << "," << stats.leaves << ","
						<< stats.pruned_by_bound << "," << stats.pruned_by_memo_exact << ","
						<< stats.pruned_by_memo_lb << "," << stats.memo_hits << "," << stats.memo_misses << ","
						<< stats.memo_inserts << "," << stats.memo_updates << "," << stats.memo_evictions << ","
						<< stats.memo_peak_size << "," << stats.memo_final_size << ","
						<< status << "," << csv_escape(error_text) << "\n";
					csv.flush();

					std::cout << "[bench] n=" << n << " instance=" << instance_idx
						<< " due=" << to_string(base.due_model)
						<< " memo=" << to_string(mode) << " clean=" << to_string(cfg.memo_cleaning)
						<< " branch=" << to_string(branch)
						<< " status=" << status << " cost=" << cost
						<< " time_ms=" << stats.elapsed_ms << "\n";
				}
			}
		}
	}

	std::signal(SIGINT, SIG_DFL);

	if (g_stop_requested) {
		std::cout << "[bench] Interrupted. Partial CSV results are saved.\n";
	}
	return true;
}

bool run_self_check(const cli_options& base) {
	std::cout << "[self-check] Running small exact checks (DFS+B&B+memo vs brute force)\n";

	instance hand;
	hand.jobs = {
		{3, 4},
		{2, 6},
		{4, 7},
		{3, 5},
		{2, 8},
	};

	std::vector<instance> tests;
	tests.push_back(hand);
	for (int k = 0; k < 8; ++k) {
		tests.push_back(generate_random_instance(8, 1, 12, 1, 35, base.config.zobrist_seed + k + 17));
	}

	const memo_mode modes[] = {memo_mode::passive, memo_mode::predictive, memo_mode::solution};
	for (memo_mode mode : modes) {
		dfs_config cfg = base.config;
		cfg.memo = mode;
		cfg.memo_capacity = std::max<std::size_t>(cfg.memo_capacity, 5000000);

		for (std::size_t i = 0; i < tests.size(); ++i) {
			const instance& inst = tests[i];
			const long long exact = brute_force_optimal(inst, nullptr);

			dfs_solver solver(cfg);
			const solve_result res = solver.solve(inst);
			if (res.best.cost != static_cast<schedule_cost_t>(exact)) {
				std::cout << "[self-check] FAILED mode=" << to_string(mode) << " case=" << i
					<< " expected=" << exact << " got=" << res.best.cost << "\n";
				return false;
			}
		}
	}

	std::cout << "[self-check] PASSED\n";
	return true;
}

void print_solution(const instance& inst, const solve_result& result, const dfs_config& cfg, const cli_options& opts) {
	std::cout << "n: " << inst.jobs.size() << "\n";
	std::cout << "best_sum_tardiness: " << result.best.cost << "\n";
	std::cout << "due_model: " << to_string(opts.due_model) << "\n";
	if (opts.due_model == due_date_model::potts) {
		std::cout << "due_range: " << opts.due_range << "\n";
		std::cout << "due_tardiness: " << opts.due_tardiness << "\n";
	}
	std::cout << "memo_mode: " << to_string(cfg.memo) << "\n";
	std::cout << "memo_cleaning: " << to_string(cfg.memo_cleaning) << "\n";
	std::cout << "branch_rule: " << to_string(cfg.branching) << "\n";
	std::cout << "memo_capacity: " << cfg.memo_capacity << "\n";
	std::cout << "zobrist_seed: " << cfg.zobrist_seed << "\n";
	std::cout << "time_ms: " << result.stats.elapsed_ms << "\n";
	std::cout << "nodes: " << result.stats.nodes << "\n";
	std::cout << "leaves: " << result.stats.leaves << "\n";
	std::cout << "pruned_by_bound: " << result.stats.pruned_by_bound << "\n";
	std::cout << "pruned_by_memo_exact: " << result.stats.pruned_by_memo_exact << "\n";
	std::cout << "pruned_by_memo_lb: " << result.stats.pruned_by_memo_lb << "\n";
	std::cout << "memo_hits: " << result.stats.memo_hits << "\n";
	std::cout << "memo_misses: " << result.stats.memo_misses << "\n";
	std::cout << "memo_inserts: " << result.stats.memo_inserts << "\n";
	std::cout << "memo_updates: " << result.stats.memo_updates << "\n";
	std::cout << "memo_evictions: " << result.stats.memo_evictions << "\n";
	std::cout << "memo_peak_size: " << result.stats.memo_peak_size << "\n";
	std::cout << "memo_final_size: " << result.stats.memo_final_size << "\n";

	if (!result.best.order.empty()) {
		std::cout << "best_order(job_indices):";
		for (int j : result.best.order) {
			std::cout << " " << j;
		}
		std::cout << "\n";
	}
}
} // namespace

int main(int argc, char** argv) {
	cli_options opts;
	opts.n = 30;
	opts.n_from = opts.n;
	opts.n_to = opts.n;
	opts.config.memo = memo_mode::predictive;
	opts.config.memo_capacity = 50000000;
	opts.config.branching = branch_rule::slack;
	opts.config.zobrist_seed = 1;
	std::string error;
	if (!parse_cli(argc, argv, opts, error)) {
		std::cerr << "Argument error: " << error << "\n\n";
		print_usage();
		return 1;
	}
	if (opts.show_help) {
		print_usage();
		return 0;
	}

	try {
		if (opts.self_check) {
			return run_self_check(opts) ? 0 : 2;
		}
		if (opts.benchmark) {
			std::string bench_error;
			if (!run_benchmark_csv(opts, bench_error)) {
				std::cerr << "Benchmark error: " << bench_error << "\n";
				return 2;
			}
			return 0;
		}

		instance inst;
		if (!opts.file_path.empty()) {
			std::ifstream file(opts.file_path);
			if (!file) {
				std::cerr << "Failed to open file: " << opts.file_path << "\n";
				return 2;
			}
			std::string parse_error;
			if (!parse_instance(file, inst, &parse_error)) {
				std::cerr << "Failed to load instance: " << parse_error << "\n";
				return 2;
			}
		}
		else {
			inst = generate_instance_for_run(opts, opts.n, opts.config.zobrist_seed, opts.d_max);
		}

		if (opts.config.memo == memo_mode::none) {
			opts.config.memo_capacity = 0;
		}
		else if (!opts.memo_cap_explicit) {
			opts.config.memo_capacity = recommended_memo_capacity(inst.jobs.size(), opts.memory_budget_mb);
		}

		dfs_solver solver(opts.config);
		const solve_result result = solver.solve(inst);
		print_solution(inst, result, opts.config, opts);
		return 0;
	}
	catch (const std::exception& ex) {
		std::cerr << "Error: " << ex.what() << "\n";
		return 3;
	}
}
