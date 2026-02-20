#pragma once

#include <cstdint>

#include "core.h"

instance generate_random_instance(int n, int p_min, int p_max, int d_min, int d_max, std::uint64_t seed);
instance generate_potts_instance(int n, int p_min, int p_max, double due_range, double tardiness_factor,
	std::uint64_t seed);
