#pragma once
#include <cstdint>
#include <vector>
#include <istream>
#include <string>

using processing_time_t = std::uint16_t;	
using due_date_t = std::uint16_t;
using schedule_time_t = std::uint64_t;	// C
using schedule_cost_t = std::uint64_t;	// Sum-Tj

struct job {
	processing_time_t p = 0;
	due_date_t d = 0;
};

struct instance {
	std::vector<job> jobs;
};

struct schedule {
	std::vector<int> order;
	schedule_cost_t cost = 0;
};

inline schedule_cost_t tardiness(schedule_time_t completion_time, due_date_t due_date) {
	return (completion_time > due_date) ? (completion_time - due_date) : 0;
}

inline schedule_cost_t evaluate_sum_tardiness(const instance& inst, const std::vector<int>& order) {
	schedule_time_t t = 0;
	schedule_cost_t sum = 0;
	for (int j : order) {
		t += inst.jobs[static_cast<std::size_t>(j)].p;
		sum += tardiness(t, inst.jobs[static_cast<std::size_t>(j)].d);
	}
	return sum;
}

bool parse_instance(std::istream& in, instance& out, std::string* error = nullptr);
bool validate_instance(const instance& inst, std::string* error = nullptr);
