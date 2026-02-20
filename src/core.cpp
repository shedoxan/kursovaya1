#include "core.h"
#include <limits>
#include <sstream>


static bool parse_job_line(const std::string& raw_line, int line_no, job& out, std::string* error) {
	std::istringstream ss(raw_line);
	long long p_value = 0;
	long long d_value = 0;
	std::string extra;
	if (!(ss >> p_value >> d_value)) {
		if (error != nullptr) {
			*error = "Invalid job format at line " + std::to_string(line_no) + ". Expected 'p d'.";
		}
		return false;
	}
	if (ss >> extra) {
		if (error != nullptr) {
			*error = "Invalid job format at line " + std::to_string(line_no) + ". Expected exactly two integers: 'p d'.";
		}
		return false;
	}

	if (p_value < std::numeric_limits<processing_time_t>::min() ||
		p_value > std::numeric_limits<processing_time_t>::max() ||
		d_value < std::numeric_limits<due_date_t>::min() ||
		d_value > std::numeric_limits<due_date_t>::max()) {
		if (error != nullptr) {
			*error = "Job value is out of supported 16-bit unsigned range at line " + std::to_string(line_no) + ".";
		}
		return false;
	}

	out.p = static_cast<processing_time_t>(p_value);
	out.d = static_cast<due_date_t>(d_value);
	return true;
}

bool parse_instance(std::istream& in, instance& out, std::string* error) {
	out.jobs.clear();

	std::string line;
	int line_no = 0;
	int n = -1;

	while (std::getline(in, line)) {
		++line_no;
		if (line.empty()) {
			continue;
		}

		std::istringstream ss(line);
		long long n_value = -1;
		std::string extra;
		if (!(ss >> n_value) || (ss >> extra) || n_value < 0 || n_value > std::numeric_limits<int>::max()) {
			if (error != nullptr) {
				*error = "Expected non-negative integer n on line " + std::to_string(line_no) + ".";
			}
			return false;
		}
		n = static_cast<int>(n_value);
		break;
	}

	if (n < 0) {
		if (error != nullptr) {
			*error = "Input does not contain the number of jobs.";
		}
		return false;
	}

	out.jobs.reserve(static_cast<std::size_t>(n));
	while (static_cast<int>(out.jobs.size()) < n && std::getline(in, line)) {
		++line_no;

		if (line.empty()) {
			continue;
		}

		job current{};
		if (!parse_job_line(line, line_no, current, error)) {
			return false;
		}
		out.jobs.push_back(current);
	}

	if (static_cast<int>(out.jobs.size()) != n) {
		if (error != nullptr) {
			*error = "Expected " + std::to_string(n) + " jobs but parsed " +
				std::to_string(out.jobs.size()) + ".";
		}
		return false;
	}

	return validate_instance(out, error);
}

bool validate_instance(const instance& inst, std::string* error) {
	long long total_processing_time = 0;
	for (std::size_t i = 0; i < inst.jobs.size(); ++i) {
		const job& j = inst.jobs[i];
		if (j.p <= 0) {
			if (error != nullptr) {
				*error = "Job at index " + std::to_string(i) + " has non-positive processing time.";
			}
			return false;
		}
		total_processing_time += j.p;
		if (total_processing_time > std::numeric_limits<int>::max()) {
			if (error != nullptr) {
				*error = "Total processing time exceeds supported limit for the current solver implementation.";
			}
			return false;
		}
	}
	return true;
}