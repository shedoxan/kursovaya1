#include "generator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

instance generate_random_instance(int n, int p_min, int p_max, int d_min, int d_max, std::uint64_t seed) {
	if (n < 0) {
		throw std::invalid_argument("n must be non-negative.");
	}
	if (p_min <= 0 || p_max < p_min) {
		throw std::invalid_argument("Invalid processing-time range.");
	}
	if (p_max > static_cast<int>(std::numeric_limits<processing_time_t>::max())) {
		throw std::invalid_argument("Processing-time range exceeds supported type limits.");
	}
	if (d_min < 0) {
		throw std::invalid_argument("Due-date range cannot be negative.");
	}
	if (d_max < d_min) {
		throw std::invalid_argument("Invalid due-date range.");
	}
	if (d_max > static_cast<int>(std::numeric_limits<due_date_t>::max())) {
		throw std::invalid_argument("Due-date range exceeds supported type limits.");
	}

	std::mt19937_64 rng(seed);
	std::uniform_int_distribution<int> p_dist(p_min, p_max);
	std::uniform_int_distribution<int> d_dist(d_min, d_max);

	instance inst;
	inst.jobs.reserve(static_cast<std::size_t>(n));
	for (int i = 0; i < n; ++i) {
		job j{};
		j.p = static_cast<processing_time_t>(p_dist(rng));
		j.d = static_cast<due_date_t>(d_dist(rng));
		inst.jobs.push_back(j);
	}

	return inst;
}

instance generate_potts_instance(int n, int p_min, int p_max, double due_range, double tardiness_factor,
	std::uint64_t seed) {
	if (n < 0) {
		throw std::invalid_argument("n must be non-negative.");
	}
	if (p_min <= 0 || p_max < p_min) {
		throw std::invalid_argument("Invalid processing-time range.");
	}
	if (p_max > static_cast<int>(std::numeric_limits<processing_time_t>::max())) {
		throw std::invalid_argument("Processing-time range exceeds supported type limits.");
	}
	if (!std::isfinite(due_range) || !std::isfinite(tardiness_factor)) {
		throw std::invalid_argument("Invalid due-date model parameters.");
	}
	if (due_range <= 0.0) {
		throw std::invalid_argument("--due-range must be positive.");
	}
	if (tardiness_factor < 0.0 || tardiness_factor > 1.5) {
		throw std::invalid_argument("--due-tardiness must be in [0, 1.5].");
	}

	std::mt19937_64 rng(seed);
	std::uniform_int_distribution<int> p_dist(p_min, p_max);

	instance inst;
	inst.jobs.reserve(static_cast<std::size_t>(n));
	long long total_p = 0;
	for (int i = 0; i < n; ++i) {
		job j{};
		j.p = static_cast<processing_time_t>(p_dist(rng));
		total_p += j.p;
		inst.jobs.push_back(j);
	}

	const double u = 1.0 - tardiness_factor - due_range / 2.0;
	const double v = 1.0 - tardiness_factor + due_range / 2.0;
	const double low_raw = static_cast<double>(total_p) * std::min(u, v);
	const double high_raw = static_cast<double>(total_p) * std::max(u, v);
	const int low_due = static_cast<int>(std::floor(low_raw));
	const int high_due = static_cast<int>(std::ceil(high_raw));
	const int clamped_low = std::max(low_due, 0);
	const int clamped_high = std::max(clamped_low, high_due);
	const int max_due = static_cast<int>(std::numeric_limits<due_date_t>::max());

	std::uniform_int_distribution<int> d_dist(std::min(clamped_low, max_due), std::min(clamped_high, max_due));
	for (job& j : inst.jobs) {
		j.d = static_cast<due_date_t>(d_dist(rng));
	}

	return inst;
}
