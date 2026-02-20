#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <list>
#include <map>
#include <unordered_map>
#include <vector>

#include "solver.h"

struct memo_lookup_result {
	bool found = false;
	bool has_exact = false;
	long long exact = 0;
	long long lower_bound = 0;
	int best_job = -1;
};

struct memo_table_stats {
	std::uint64_t hits = 0;
	std::uint64_t misses = 0;
	std::uint64_t inserts = 0;
	std::uint64_t updates = 0;
	std::uint64_t evictions = 0;
	std::size_t peak_size = 0;
	std::size_t final_size = 0;
};

class memo_table {
public:
	explicit memo_table(std::size_t capacity = 0,
		memo_cleaning_policy cleaning = memo_cleaning_policy::lru)
		: capacity_(capacity), cleaning_(cleaning) {}

	void clear() {
		entries_.clear();
		buckets_.clear();
		lufo_order_.clear();
		lufo_by_id_.clear();
		next_id_ = 1;
		stats_ = {};
	}

	void set_capacity(std::size_t capacity, bool count_stats = true) {
		capacity_ = capacity;
		while (entries_.size() > capacity_) {
			evict_one(count_stats);
		}
		stats_.final_size = entries_.size();
	}

	void set_cleaning_policy(memo_cleaning_policy policy, bool count_stats = true) {
		if (cleaning_ == policy) {
			return;
		}
		cleaning_ = policy;
		if (cleaning_ == memo_cleaning_policy::lufo) {
			build_lufo_index();
		}
		else {
			clear_lufo_index();
		}
		while (entries_.size() > capacity_) {
			evict_one(count_stats);
		}
		stats_.final_size = entries_.size();
	}

	[[nodiscard]] std::size_t size() const {
		return entries_.size();
	}

	[[nodiscard]] std::size_t capacity() const {
		return capacity_;
	}

	[[nodiscard]] memo_lookup_result lookup(const std::vector<std::uint64_t>& bits, int time,
		std::uint64_t subset_hash, bool count_stats = true) {
		const std::uint64_t key_hash = make_state_hash(subset_hash, time);
		list_it* it = find_entry(bits, time, key_hash);
		if (it == nullptr) {
			if (count_stats) {
				++stats_.misses;
			}
			return {};
		}

		if (count_stats) {
			++stats_.hits;
		}

		note_use(*it);

		memo_lookup_result out;
		out.found = true;
		out.has_exact = (*it)->has_exact;
		out.exact = (*it)->exact;
		out.lower_bound = (*it)->lower_bound;
		out.best_job = (*it)->best_job;
		return out;
	}

	void store_lower_bound(const std::vector<std::uint64_t>& bits, int time,
		std::uint64_t subset_hash, long long lb, bool count_stats = true) {
		if (capacity_ == 0) {
			return;
		}

		const std::uint64_t key_hash = make_state_hash(subset_hash, time);
		list_it* it = find_entry(bits, time, key_hash);
		if (it != nullptr) {
			bool changed = false;
			if (lb > (*it)->lower_bound) {
				(*it)->lower_bound = lb;
				changed = true;
			}
			note_touch(*it);
			if (count_stats && changed) {
				++stats_.updates;
			}
			return;
		}

		ensure_capacity_for_insert(count_stats);

		entry e;
		e.bits = bits;
		e.time = time;
		e.state_hash = key_hash;
		e.lower_bound = lb;
		e.has_exact = false;
		e.exact = 0;
		e.best_job = -1;
		e.use_count = 0;
		e.id = next_id_++;
		e.has_lufo_pos = false;

		entries_.push_front(std::move(e));
		const list_it inserted = entries_.begin();
		buckets_[key_hash].push_back(inserted);
		if (cleaning_ == memo_cleaning_policy::lufo) {
			register_lufo_entry(inserted);
		}

		if (count_stats) {
			++stats_.inserts;
		}
		stats_.peak_size = std::max(stats_.peak_size, entries_.size());
		stats_.final_size = entries_.size();
	}

	void store_exact(const std::vector<std::uint64_t>& bits, int time,
		std::uint64_t subset_hash, long long exact, int best_job, bool count_stats = true) {
		if (capacity_ == 0) {
			return;
		}

		const std::uint64_t key_hash = make_state_hash(subset_hash, time);
		list_it* it = find_entry(bits, time, key_hash);
		if (it != nullptr) {
			bool changed = false;
			if (!(*it)->has_exact || exact < (*it)->exact) {
				(*it)->has_exact = true;
				(*it)->exact = exact;
				(*it)->best_job = best_job;
				changed = true;
			}
			if (exact > (*it)->lower_bound) {
				(*it)->lower_bound = exact;
				changed = true;
			}
			note_touch(*it);
			if (count_stats && changed) {
				++stats_.updates;
			}
			return;
		}

		ensure_capacity_for_insert(count_stats);

		entry e;
		e.bits = bits;
		e.time = time;
		e.state_hash = key_hash;
		e.lower_bound = exact;
		e.has_exact = true;
		e.exact = exact;
		e.best_job = best_job;
		e.use_count = 0;
		e.id = next_id_++;
		e.has_lufo_pos = false;

		entries_.push_front(std::move(e));
		const list_it inserted = entries_.begin();
		buckets_[key_hash].push_back(inserted);
		if (cleaning_ == memo_cleaning_policy::lufo) {
			register_lufo_entry(inserted);
		}

		if (count_stats) {
			++stats_.inserts;
		}
		stats_.peak_size = std::max(stats_.peak_size, entries_.size());
		stats_.final_size = entries_.size();
	}

	[[nodiscard]] memo_table_stats stats() const {
		memo_table_stats out = stats_;
		out.final_size = entries_.size();
		return out;
	}

private:
	using lufo_it = std::multimap<std::uint64_t, std::uint64_t>::iterator;

	struct entry {
		std::vector<std::uint64_t> bits;
		int time = 0;
		std::uint64_t state_hash = 0;
		long long lower_bound = 0;
		bool has_exact = false;
		long long exact = 0;
		int best_job = -1;
		std::uint64_t use_count = 0;
		std::uint64_t id = 0;
		lufo_it lufo_pos{};
		bool has_lufo_pos = false;
	};

	static std::uint64_t mix64(std::uint64_t x) {
		x += 0x9e3779b97f4a7c15ULL;
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
		x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
		return x ^ (x >> 31);
	}

	static std::uint64_t make_state_hash(std::uint64_t subset_hash, int time) {
		return subset_hash ^ mix64(static_cast<std::uint64_t>(static_cast<std::uint32_t>(time)));
	}

	using list_it = std::list<entry>::iterator;

	list_it* find_entry(const std::vector<std::uint64_t>& bits, int time, std::uint64_t key_hash) {
		auto bucket_it = buckets_.find(key_hash);
		if (bucket_it == buckets_.end()) {
			return nullptr;
		}

		auto& vec = bucket_it->second;
		for (auto it = vec.begin(); it != vec.end(); ++it) {
			if ((**it).time == time && (**it).bits == bits) {
				return &(*it);
			}
		}
		return nullptr;
	}

	void note_use(list_it it) {
		if (cleaning_ == memo_cleaning_policy::lru) {
			entries_.splice(entries_.begin(), entries_, it);
			return;
		}

		if (!it->has_lufo_pos) {
			register_lufo_entry(it);
		}
		else {
			lufo_order_.erase(it->lufo_pos);
		}
		++(it->use_count);
		it->lufo_pos = lufo_order_.insert({it->use_count, it->id});
		it->has_lufo_pos = true;
		lufo_by_id_[it->id] = it;
	}

	void note_touch(list_it it) {
		if (cleaning_ == memo_cleaning_policy::lru) {
			entries_.splice(entries_.begin(), entries_, it);
		}
	}

	void ensure_capacity_for_insert(bool count_stats) {
		while (entries_.size() >= capacity_ && !entries_.empty()) {
			evict_one(count_stats);
		}
	}

	void evict_one(bool count_stats) {
		if (entries_.empty()) {
			return;
		}
		if (cleaning_ == memo_cleaning_policy::lufo) {
			evict_lufo_min_batch(count_stats);
			return;
		}
		evict_lru_tail(count_stats);
	}

	void evict_lru_tail(bool count_stats) {
		if (entries_.empty()) {
			return;
		}
		erase_entry(std::prev(entries_.end()), count_stats);
	}

	void evict_lufo_min_batch(bool count_stats) {
		if (lufo_order_.empty()) {
			evict_lru_tail(count_stats);
			return;
		}

		const std::uint64_t min_use = lufo_order_.begin()->first;
		bool removed_any = false;
		while (!lufo_order_.empty() && lufo_order_.begin()->first == min_use) {
			const std::uint64_t id = lufo_order_.begin()->second;
			auto id_it = lufo_by_id_.find(id);
			if (id_it == lufo_by_id_.end()) {
				lufo_order_.erase(lufo_order_.begin());
				continue;
			}
			erase_entry(id_it->second, count_stats);
			removed_any = true;
		}

		if (!removed_any) {
			evict_lru_tail(count_stats);
		}
	}

	void erase_entry(list_it victim, bool count_stats) {
		const std::uint64_t key_hash = victim->state_hash;

		auto bucket_it = buckets_.find(key_hash);
		if (bucket_it != buckets_.end()) {
			auto& vec = bucket_it->second;
			vec.erase(std::remove(vec.begin(), vec.end(), victim), vec.end());
			if (vec.empty()) {
				buckets_.erase(bucket_it);
			}
		}

		if (victim->has_lufo_pos) {
			lufo_order_.erase(victim->lufo_pos);
			lufo_by_id_.erase(victim->id);
		}

		entries_.erase(victim);
		if (count_stats) {
			++stats_.evictions;
		}
		stats_.final_size = entries_.size();
	}

	void register_lufo_entry(list_it it) {
		it->lufo_pos = lufo_order_.insert({it->use_count, it->id});
		it->has_lufo_pos = true;
		lufo_by_id_[it->id] = it;
	}

	void build_lufo_index() {
		clear_lufo_index();
		for (auto it = entries_.begin(); it != entries_.end(); ++it) {
			register_lufo_entry(it);
		}
	}

	void clear_lufo_index() {
		lufo_order_.clear();
		lufo_by_id_.clear();
		for (auto it = entries_.begin(); it != entries_.end(); ++it) {
			it->has_lufo_pos = false;
		}
	}

	std::size_t capacity_ = 0;
	memo_cleaning_policy cleaning_ = memo_cleaning_policy::lru;
	std::list<entry> entries_;
	std::unordered_map<std::uint64_t, std::vector<list_it>> buckets_;
	std::multimap<std::uint64_t, std::uint64_t> lufo_order_;
	std::unordered_map<std::uint64_t, list_it> lufo_by_id_;
	std::uint64_t next_id_ = 1;
	memo_table_stats stats_{};
};
