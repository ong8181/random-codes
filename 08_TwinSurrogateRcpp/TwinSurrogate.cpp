#include <Rcpp.h>
using namespace Rcpp;

#include <array>

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericMatrix cpp_dist_euclidean(const NumericMatrix& inMat)
{
	int	tNRows = inMat.rows();
	int	tNCols = inMat.cols();
	NumericMatrix	outMat = NumericMatrix(tNRows, tNRows);
	for (int i = 0; i < tNRows - 1; ++i) {
		const auto&	tRow1 = inMat.row(i);
		for (int j = i + 1; j < tNRows; ++j) {
			const auto&	tRow2 = inMat.row(j);
			// outMat(i, j) = outMat(j, i) = sqrt(sum(pow(tRow1 - tRow2, 2)));
			double	s = 0;
			for (int c = 0; c < tNCols; ++c) {
				double d = tRow1[c] - tRow2[c];
				s += d * d;
			}
			outMat(i, j) = outMat(j, i) = std::sqrt(s);
		}
	}
	return outMat;
}

// [[Rcpp::export]]
NumericMatrix cpp_dist_maxnorm(const NumericMatrix& inMat)
{
	int	tNRows = inMat.rows();
	int	tNCols = inMat.cols();
	NumericMatrix	outMat = NumericMatrix(tNRows, tNRows);
	for (int i = 0; i < tNRows - 1; ++i) {
		const auto&	tRow1 = inMat.row(i);
		for (int j = i + 1; j < tNRows; ++j) {
			const auto&	tRow2 = inMat.row(j);
			// outMat(i, j) = outMat(j, i) = max(abs(tRow1 - tRow2));
			double m = 0;
			for (int c = 0; c < tNCols; ++c) {
				double d = tRow1[c] - tRow2[c];
				if (d < 0) {
					d = -d;
				}
				if (d > m) {
					m = d;
				}
			}
			outMat(i, j) = outMat(j, i) = m;
		}
	}
	return outMat;
}

// [[Rcpp::export]]
double cpp_quantile(const NumericVector& inVec, const double q)
{
	int n = inVec.length();
	double index = (double)(n - 1) * q;
	int i_lo = std::floor(index);
	int i_hi = std::ceil(index);
	NumericVector tSorted = clone(inVec);
	std::sort(tSorted.begin(), tSorted.end());
	double q_lo = tSorted[i_lo];
	double q_hi = tSorted[i_hi];
	if (q_lo != q_hi) {
		double h = index - i_lo;
		return q_lo * (1.0 - h) + q_hi * h;
	} else {
		return q_lo;
	}
}

// [[Rcpp::export]]
IntegerMatrix cpp_greater(const NumericMatrix& inMat, const double inThreshold)
{
	// return IntegerMatrix(inMat.nrow(), inMat.ncol(),
	//						ifelse(inMat > inThreshold, 1, 0).begin());
	IntegerMatrix m(inMat.nrow(), inMat.ncol());
	auto iter = m.begin();
	for (const auto srciter : inMat) {
		if (srciter > inThreshold) {
			*iter = 1;
		// } else {
		// 	*iter = 0;
		}
		++iter;
	}
	return m;
}

// [[Rcpp::export]]
IntegerMatrix cpp_distOneZero(const NumericMatrix& inMat, const std::string& inMethod, double s = 0.875)
{
	BEGIN_RCPP

	NumericMatrix dist_mat;
	if (inMethod == "norm") {
		dist_mat = cpp_dist_maxnorm(inMat);
	} else { //}if (inMethod == "euclid")
		dist_mat = cpp_dist_euclidean(inMat);
	}
	double d_threshold = cpp_quantile(dist_mat, s);
	return cpp_greater(dist_mat, d_threshold);

	END_RCPP	
}

int PointNext(int x, const std::vector<std::pair<int, int> >& twins, int d_ncol)
{
	int nex;
	if (x == -1) {
		nex = -1;
	} else {
		std::vector<int> twins_x;
		for (const auto& t : twins) {
			if (t.first == x) {
				twins_x.push_back(t.second);
			}
		}
		if (twins_x.size() > 0) {
			nex = twins_x[std::floor(runif(1)[0] * twins_x.size())] + 1;
		} else {
			nex = x + 1;
		}
	}
	if (nex >= d_ncol) {
		nex = -1;
	}
	return nex;
}

// [[Rcpp::export]]
NumericMatrix cpp_twinSurrogate(const NumericMatrix& original_e,
                                int dim,
                                int num_iter,
                                double s = 0.875,
                                const std::string& surrogate_option = "random",
                                const std::string& initial_point = "same_season",
                                const std::string& distance_method = "norm",
                                int point_per_year = 24,
                                const std::string& s_update = "on",
                                int n_twin_threshold = 10,
                                bool output_message = false)
{
	BEGIN_RCPP
	std::vector<double> s_seq = {0.875, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93,
								 0.94, 0.95, 0.87, 0.86, 0.85, 0.84, 0.83,
								 0.82, 0.81, 0.80};
	s_seq[0] = s;

	IntegerMatrix d;
	std::vector<std::pair<int, int>> twins;
	bool s_update_is_off = (s_update == "off");
	bool surrogate_option_is_phase_lock = (surrogate_option == "phase_lock");
	int num_twins = 0;
	for (auto s : s_seq) {
		d = cpp_distOneZero(original_e, distance_method, s);
		
		twins.clear();
		size_t d_ncol = d.ncol();
		size_t d_nrow = d.nrow();
		for (int i = 0; i < d_ncol; ++i) {
			for (int j = 0; j < d_ncol; ++j) {
				int k = 0;
				for (; k < d_nrow; ++k) {
					if (d(k, i) != d(k, j)) {
						break;
					}
				}
				if (k == d_nrow) {
					if (surrogate_option_is_phase_lock) {
						if ((i - j) % point_per_year == 0) {
							twins.emplace_back(i, j);
						}
					} else {
						twins.emplace_back(i, j);
					}
				}
			}
		}
		num_twins = twins.size() - d_nrow;
		if (num_twins >= n_twin_threshold) {
			break;
		}
		if (s_update_is_off) {
			break;
		}
	}

	size_t d_ncol = d.ncol();
	size_t d_nrow = d.nrow();
	if (output_message) {
		double prop_black = (double)sum(d) / (double)(d_ncol * d_nrow);
		std::cout << "Proportion of 1:" << prop_black << std::endl;
		std::cout << "Number of twins:" << num_twins << std::endl;
	}

	NumericMatrix surrogate(dim - 1 + d_ncol, num_iter);
	int col = 0;
	int avoid_infinite_loop = 0;
	while (1) {
		int surr0 = 0;
		if (surrogate_option == "random") {
			// Select random initial points
			surr0 = std::floor(runif(1)[0] * (d_ncol - 1));
		} else if (surrogate_option == "phase_lock") {
			if (initial_point == "same_season") {
				// Select the point of the same season
				surr0 = point_per_year + point_per_year *
							std::floor(runif(1)[0] *
							(int)((d_ncol - 2) / point_per_year));
			} else if (initial_point == "twins") {
				// Select twins of the original initial point as the surrogate initial point
				std::vector<int> twins_zero;
				for (const auto& t : twins) {
					if (t.first == 0) {
						twins_zero.push_back(t.second);
					}
				}
				surr0 = twins_zero[std::floor(runif(1)[0] * twins_zero.size())];
			}
		} else {
			std::cerr << "Warning: specify the correct option!" << std::endl;
		}
		
		// Search next point
		std::vector<int> surr(d_ncol);
		surr[0] = surr0;
		for (int j = 0; j < d_ncol - 1; ++j) {
			int nex;
			while (1) {
				nex = PointNext(surr[j], twins, d_ncol);
				if (surrogate_option == "phase_lock" &&
					initial_point == "same_season") {
					if (nex != j + 1) {
						break;
					}
				} else {
					break;
				}
			}
			surr[j + 1] = nex;
		}
		// Save the surrogate if it reach to the length of the origial data
		// Not save if the surrogate is short
		if (surr.back() != -1) {
			NumericMatrix::Column tColumn = surrogate.column(col);
			NumericMatrix::ConstRow tRow = original_e.row(surr[0]);
			int i = 0;
			for (; i < dim - 1; ++i) {
				tColumn[i] = tRow[i];
			}
			for (int r : surr) {
				// std::cerr << ".." << original_e(r, dim - 1) << std::endl;
				tColumn[i++] = original_e(r, dim - 1);
			}
			// std::cerr << i << ',' << dim - 1 + d_ncol << std::endl;
			++col;
		}

		// Complete the surrogate generation if length(surrogate) == num.iter
		if (col >= num_iter) {
			break;
		}

		// Avoid infinite loop
		// End cycles if surrogates cannot be generated during > 30*num.iter trials
		if (++avoid_infinite_loop > 30 * num_iter) {
			break;
		}
	}

	if (avoid_infinite_loop <= 30 * num_iter) {
		return surrogate;
	} else {
		NumericVector v = NumericVector((dim - 1 + d_ncol) * num_iter, R_NaN);
		v.attr("dim") = Dimension(dim - 1 + d_ncol, num_iter);
		return as<NumericMatrix>(v);
	}
	END_RCPP
}

// [[Rcpp::export]]
NumericMatrix cpp_nanMatrix(int nrow, int ncol)
{
	NumericVector v = NumericVector(nrow * ncol, R_NaN);
	v.attr("dim") = Dimension(nrow, ncol);
	return as<NumericMatrix>(v);
}

/*** R

TwinSurrogateRcpp <- function(original,
							  dim,
							  num.iter,
							  tau = 1,
							  s   = 0.875,
							  surrogate.option = c("random", "phase_lock"),
							  initial.point    = c("same_season", "twins"),
							  distance.method  = c("norm", "euclid"),
							  point.per.year   = 24,
							  s.update         = c("on", "off"),
							  n.twin.threshold = 10,
							  output.message = F) {
	surrogate.option = match.arg(surrogate.option)
	initial.point = match.arg(initial.point)
	distance.method = match.arg(distance.method)
	s.update = match.arg(s.update)
	
	# generate time-lag embedding matrix
	if (dim >  1) {
		original_e <- tseriesChaos::embedd(original, dim, tau)
	}
	if (dim == 1) {
		original_e <- as.matrix(original)
	}
	if (dim <  1) {
		cat("Warning: Embedding dimension should be >= 1")
	}
	surrogate.one.col <- cpp_twinSurrogate(
		original_e,
		dim,
		num.iter,
		s,
		surrogate.option,
		initial.point,
		distance.method,
		point.per.year,
		s.update,
		n.twin.threshold,
		output.message
	)
	return(as.data.frame(surrogate.one.col))
}

*/
