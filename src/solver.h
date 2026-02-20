#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "core.h"
#include "stats.h"

enum class memo_mode {
	none = 0,
	passive = 1,
	predictive = 2,
	solution = 3
};

enum class memo_cleaning_policy {
	lru = 0,
	lufo = 1
};

enum class branch_rule {
	slack = 0,
	edd = 1,
	spt = 2,
	decomposition = 3
};

struct dfs_config {
	memo_mode memo = memo_mode::predictive;
	branch_rule branching = branch_rule::slack;
	memo_cleaning_policy memo_cleaning = memo_cleaning_policy::lru;
	std::size_t memo_capacity = 200000;
	std::uint64_t zobrist_seed = 1;
	bool reconstruct_order = true;
};

const char* to_string(memo_mode mode);
const char* to_string(memo_cleaning_policy policy);
const char* to_string(branch_rule rule);
bool parse_memo_mode(const std::string& text, memo_mode& out);
bool parse_memo_cleaning_policy(const std::string& text, memo_cleaning_policy& out);
bool parse_branch_rule(const std::string& text, branch_rule& out);

struct solve_result {
	schedule best;
	solver_stats stats;
};

class solver {
public:
	virtual ~solver() = default;
	virtual solve_result solve(const instance& inst) = 0;
};
