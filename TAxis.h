#ifndef TAXIS_H_
#define TAXIS_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "utilities.h"

template<class T>
class TAxis {

private:
	T min;
	T max;
	T reference_value;
	size_t size;
	size_t idx;
	vector<T> axis;

public:
	TAxis() :
		min(0), max(0), reference_value(0), size(0), idx(0) {
	}

	virtual ~TAxis() {
		axis.clear();
	}

	/* Accessor / Mutator */
	T &get(size_t i) {
		return axis[i];
	}

	T &get(const size_t& i) {
		return axis[i];
	}

	/* Accessor */
	const T &get(size_t i) const {
		return axis[i];
	}

	const T &get(const size_t& i) const {
		return axis[i];
	}

	T at(size_t i) {
		return axis[i];
	}

	void build_log_axis(const T& min, const T& max, const size_t& size) {
		this->min = min;
		this->max = max;
		this->size = size;

		const T delta_log = exp(log(max / min) / (size - 1));

		axis.resize(size, 0.);

		for (size_t i = 0; i < size; ++i) {
			T value = exp(log(min) + (T) i * log(delta_log));
			axis.at(i) = value;
		}
	}

	void build_lin_axis(const T& min, const T& max, const size_t& size) {
		this->min = min;
		this->max = max;
		this->size = size;

		const T dx = (max - min) / (T) (size - 1);

		axis.resize(size, 0.);

		for (size_t i = 0; i < size; ++i) {
			T value = min + dx * i;
			axis.at(i) = value;
		}
	}

	void show_axis(const string& name, const T& units) const {
		cout << name << " axis in : "
				<< *min_element(axis.begin(), axis.end()) / units << " ... ";
		cout << *max_element(axis.begin(), axis.end()) / units << " with "
				<< axis.size() << " elements.";
		cout << "\n";
	}

	T get_max() const {
		return max;
	}

	T get_min() const {
		return min;
	}

	size_t get_idx() const {
		return idx;
	}

	size_t get_idx(T value) {
		if (value < 0.999 * min || value > max * 1.001) {
			cerr << "Error: Value " << value << " cannot be set in the range : "
					<< min << " ... " << max << "\n";
			exit(-1);
		} else {
			return find_nearest(axis, value);
		}
	}

	size_t get_lower_idx(T value) {
		if (value < min || value > max) {
			cerr << "Value not in the range " << __FILE__ << __LINE__ << "\n";
			exit(-1);
		}
		size_t out = 0;
		while (axis.at(out) <= value)
			out++;
		return out - 1;
	}

	T get_reference_value() const {
		return reference_value;
	}

	void set_reference_value(T referenceValue) {
		if (referenceValue < min || referenceValue > max) {
			cerr << "Error: Reference_value " << referenceValue
					<< " cannot be set!" << "\n";
			exit(-1);
		} else {
			reference_value = referenceValue;
			idx = find_nearest(axis, reference_value);
		}
	}

	size_t get_size() const {
		return size;
	}
};

#endif /* TAXIS_H_ */
