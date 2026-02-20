#pragma once

#include <cstdint>
#include <vector>

#include "memo.h"
#include "solver.h"

class dfs_solver final : public solver {
public:
	explicit dfs_solver(dfs_config config = {});
	solve_result solve(const instance& inst) override;

private:
	long long solve_state(int depth, int current_time, const memo_lookup_result* known_lookup, bool use_exact_memo,
		bool use_lb_memo, bool track_stats);
	long long lower_bound_additional(int depth, int current_time) const;
	long long heuristic_upper_bound_edd(int depth, int current_time, int* first_job) const;
	std::vector<int> build_branch_indices(int depth, int current_time) const;
	std::vector<int> reconstruct_order(long long optimal_cost);

	void remove_job_from_state(int job_idx);
	void restore_job_to_state(int job_idx);
	bool is_job_remaining(int job_idx) const;

	dfs_config config_;
	const instance* inst_ = nullptr;
	int n_ = 0;

	std::vector<int> perm_jobs_;
	std::vector<std::uint64_t> remaining_bits_;
	std::vector<std::uint64_t> zobrist_job_;
	std::uint64_t subset_hash_ = 0;

	memo_table memo_;
	solver_stats stats_{};
};
