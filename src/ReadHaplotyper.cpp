#include "include/TandemTwister.hpp"

// #include <mlpack/core.hpp>
#include <mlpack/methods/dbscan.hpp>
#include "mlpack/namespace_compat.hpp"
#include <cassert>


float TandemTwister::calculateFeatureMeanInCluster(const std::vector<uint16_t>& cluster, const arma::mat &region_features, uint16_t dataIndex) {
    float mean = 0;
    for (uint16_t i = 0; i < cluster.size(); ++i) {
        mean += region_features(cluster[i], dataIndex);
    }
    return mean / cluster.size();
}


std::vector<uint16_t> TandemTwister::removeOutlierReadsZScore(const arma::mat &data, float zThreshold, uint16_t columnIdx) {
    std::vector<uint16_t> removedIndices;

    float mean = arma::mean(data.col(columnIdx));
    float stdDev = arma::stddev(data.col(columnIdx));

    for (uint16_t i = 0; i < data.n_rows; ++i) {
        float zScore = (data(i, columnIdx) - mean) / stdDev;

        if (std::abs(zScore) > zThreshold) {
            removedIndices.push_back(i);
        }
    }
    return removedIndices;
}

arma::mat TandemTwister::standardizeMatrix(arma::mat matrix) {
    // Calculate mean and standard deviation for each column
    arma::rowvec columnMeans = arma::mean(matrix, 0);
    arma::rowvec columnStdDevs = arma::stddev(matrix, 0, 0);

    // Standardize the matrix
    for (arma::uword col = 0; col < matrix.n_cols; col++) {
        if (columnStdDevs[col] != 0.0) {
            matrix.col(col) -= columnMeans[col];
            matrix.col(col) /= columnStdDevs[col];
        }
        else {
            matrix.col(col) -= columnMeans[col];
        }
    }
    return matrix;
}



std::tuple<std::vector<std::vector<uint16_t>>, std::vector<uint16_t>> TandemTwister::runDBSCAN(const arma::mat &features_matrix_stdandarized,float eps, uint16_t minPts)
{
    
    const double epsilon = eps;
    uint16_t minPoints = minPts;
    
    mlpack::DBSCAN<> dbscan(epsilon, minPoints);

    arma::Row<size_t> assignments;

    arma::mat transposed_features = features_matrix_stdandarized.t();
    dbscan.Cluster(transposed_features, assignments);
    
    assert(assignments.n_elem == features_matrix_stdandarized.n_rows);
    arma::Row<size_t> uniqueLabels = arma::unique(assignments);

    std::vector<uint16_t> noise;
    std::vector<std::vector<uint16_t>> clusters;
    // put the clusters in a vector of vectors
    for (uint16_t i = 0; i < uniqueLabels.n_elem; ++i)
    {
        if (uniqueLabels(i) != SIZE_MAX)
        {
            // Collect cluster points.
            std::vector<uint16_t> clusterPoints;
            for (uint16_t j = 0; j < assignments.n_elem; ++j)
            {
                if (assignments(j) == uniqueLabels(i))
                {
                    clusterPoints.push_back(j);
                }
            }
            clusters.push_back(clusterPoints);
        }
        else 
        {
            // Collect noise points.
            for (uint16_t j = 0; j < assignments.n_elem; ++j)
            {
                if (assignments(j) == SIZE_MAX)
                {
                    noise.push_back(j);
                }
            }
        }
    }
    // sort the clusters by their size form the largest to the smallest
    std::sort(clusters.begin(), clusters.end(), [](auto  a, auto b) {
        return a.size() > b.size();
    });
    // convert the clusters to vectors of ints
    return std::make_tuple(clusters, noise);

}


void TandemTwister::cluster_reads(arma::mat &region_features,std::vector<std::tuple<std::string,std::vector<Interval>,uint8_t>> &reads_intervals,std::vector<std::vector<std::pair<std::string,std::vector<Interval>>>> & final_clusters , std::vector<std::pair<std::string,std::vector<Interval>>> & noise_cluster, uint16_t motif_length){
    std::vector<std::vector<uint16_t>> clusters;
    std::vector<uint16_t> noise;
    // std::cout << "The number of reads in the region " << region_features.n_rows << std::endl;
    float eps = 0;
    unsigned int iterations = 0;
    unsigned noise_limit = 0;
    uint16_t minPts  = 0;

    minPts = std::ceil(this->minPts_frac * region_features.n_rows);
    if (motif_length <= 4){ 
        noise_limit =  std::floor(region_features.n_rows*this->noise_limit_str); 
        eps = this->start_eps_str;
    }
    else{
        noise_limit =  std::floor(region_features.n_rows*this->noise_limit_vntr);
        eps = this->start_eps_vntr;
    }
    arma::mat standardized_features;
    // standardize the features
    standardized_features = standardizeMatrix(region_features);
    float zThreshold = 2; 
    uint16_t columnIdx = 0; // The column index of the column to remove outliers from
    // if the reads type is not CCS and the number of reads is larger than 2, remove the outliers
    std::vector<uint16_t> outliers_idx = {};
    if ( this->remove_outliers_zscore && region_features.n_rows > 2) {
        outliers_idx = removeOutlierReadsZScore(region_features, zThreshold , columnIdx);
        for (uint16_t i = 0; i < outliers_idx.size(); ++i) {
            noise_cluster.push_back(std::make_pair(std::get<0>(reads_intervals[outliers_idx[i]]), std::get<1>(reads_intervals[outliers_idx[i]])));
        }
        if (!outliers_idx.empty() && outliers_idx.size() < standardized_features.n_rows) {
                // get the indices of features that are not outliers
                for (int16_t idx = standardized_features.n_rows - 1; idx >= 0; --idx) {
                    if (std::find(outliers_idx.begin(), outliers_idx.end(), idx) != outliers_idx.end()) {
                        // remove the interval from the reads_intervals and the feature row from the features
                        standardized_features.shed_row(idx);
                        region_features.shed_row(idx);
                        reads_intervals.erase(reads_intervals.begin() + idx);
                        // add it to the noise
                        noise.push_back(idx);
                    }
                }
        }

    }
    
    // if the matrix is not empty
    if (standardized_features.n_rows <= 0) {
        final_clusters.resize(2);
        // std::cout << "The number of reads in the region is 0" << std::endl;
        return;
    }
    if (standardized_features.n_rows > 1) {
        bool allColumnsEqual = true;

        for (arma::uword col = 0; col < standardized_features.n_cols; ++col) {
            if (!arma::all(standardized_features.col(col) == standardized_features(0, col))) {
                allColumnsEqual = false;
                break;
            }
        }
        if (allColumnsEqual) {
            final_clusters.resize(2);
            // std::cout << "All features are equal, not need for clustering" << std::endl;
            for (uint16_t i = 0; i < standardized_features.n_rows; ++i) {
                auto read_name = std::get<0>(reads_intervals[i]);
                auto read_intervals = std::get<1>(reads_intervals[i]);
                final_clusters[0].push_back(std::make_pair(read_name, read_intervals));
            }
            return;
        }
    }


    std::tie(clusters, noise) = runDBSCAN(standardized_features, eps, minPts);

    while (noise.size() >= noise_limit ){
        std::tie(clusters, noise) = runDBSCAN(standardized_features, eps, minPts);
        eps += 0.1;
        ++iterations;

        if (iterations >=this->cluster_iter){
            break;
        } 
    }
    std::sort(clusters.begin(), clusters.end(), [&](auto a, auto b) {
        if (a.size() == b.size()) {
            // Compare based on the read length (first value in the column/row)
            return region_features(a[0], 0) > region_features(b[0], 0);
        }
        return a.size() > b.size();
    });
    /* VISULIZE CLUSTERS */
    spdlog::debug("The number of clusters: {}", clusters.size());
    spdlog::debug("The number of noise: {}/{}", noise.size(), standardized_features.n_rows);
    spdlog::debug("The number of reads in the region: {}", region_features.n_rows);
    spdlog::debug("The number of reads in the region after removing outliers: {}", standardized_features.n_rows);


    // prnt the the clusters with their features
    for (uint16_t i = 0; i < clusters.size(); ++i) {
        spdlog::debug("Cluster {}: ", i);
        for (uint16_t j = 0; j < clusters[i].size(); ++j) {
            std::string cluster_features;
            for (uint16_t k = 0; k < standardized_features.n_cols; ++k) {
                cluster_features += fmt::format("{} ", region_features(clusters[i][j], k));
            }
            spdlog::debug("{}", cluster_features);
        }
    }
    spdlog::debug("Noise: ");
    for (uint16_t i = 0; i < noise.size();++i){
        string noise_features;
        for (uint16_t k = 0; k < standardized_features.n_cols; ++k) {
            noise_features += fmt::format("{} ", region_features(noise[i], k));
        }
        spdlog::debug("{}", noise_features);
    }
     /* VISULIZE CLUSTERS */


    // if all the reads are noise, add them to the first cluster
    if (clusters.empty()){
        std::vector<std::pair<std::string,std::vector<Interval>>> cluster_reads;
        for (uint16_t i = 0; i< standardized_features.n_rows;++i){
            auto read_name = std::get<0>(reads_intervals[i]);
            auto read_intervals = std::get<1>(reads_intervals[i]);
            cluster_reads.push_back(std::make_pair(read_name, read_intervals));
        }
        final_clusters.resize(2);
        final_clusters[0] = cluster_reads;
        return;
    }

    for (uint16_t i = 0; i < clusters.size(); ++i) {
        std::vector<std::pair<std::string,std::vector<Interval>>> cluster_reads;
        for (uint16_t j = 0; j < clusters[i].size(); ++j) {
            auto read_name = std::get<0>(reads_intervals[clusters[i][j]]);
            auto read_intervals = std::get<1>(reads_intervals[clusters[i][j]]);
            cluster_reads.push_back(std::make_pair(read_name, read_intervals));
        }
        final_clusters.push_back(cluster_reads);
    }

    for (uint16_t i = 0; i< noise.size();++i){

        auto read_name = std::get<0>(reads_intervals[noise[i]]);
        auto read_intervals = std::get<1>(reads_intervals[noise[i]]);
        noise_cluster.push_back(std::make_pair(read_name, read_intervals));
    }

    // sort the clusters by their size form the largest to the smallest
    std::sort(final_clusters.begin(), final_clusters.end(), [](auto  a, auto b) {
        return a.size() > b.size();
    });


    if (clusters.size() > 2 && this->manual_reclustering && this->somatic == false) {
        spdlog::debug("Start manual reclustering ");
        std::vector<float> first_cluster_means(region_features.n_cols, 0.0);
        std::vector<float> second_cluster_means(region_features.n_cols, 0.0);

        for (arma::uword col = 0; col < region_features.n_cols; ++col) {
            first_cluster_means[col] = calculateFeatureMeanInCluster(clusters[0], region_features, col);
            second_cluster_means[col] = calculateFeatureMeanInCluster(clusters[1], region_features, col);
        }
        float cn1_diff = 0.0;
        float cn2_diff = 0.0;
        for (uint16_t i = 2; i < clusters.size(); ++i) {
            auto& cluster = clusters[i];

            for (uint16_t j = 0; j < cluster.size(); ++j) {
                auto& read_name = std::get<0>(reads_intervals[cluster[j]]);
                auto& read_intervals = std::get<1>(reads_intervals[cluster[j]]);
                
                float avg_diff_first = 0.0;
                float avg_diff_second = 0.0;
                
                for (arma::uword col = 0; col < region_features.n_cols; ++col) {
                    float feature_value = region_features(cluster[j], col);
                    avg_diff_first += std::abs(first_cluster_means[col] - feature_value);
                    avg_diff_second += std::abs(second_cluster_means[col] - feature_value);
                    if (col == 1){
                        cn1_diff += std::abs(first_cluster_means[col] - feature_value);
                        cn2_diff += std::abs(second_cluster_means[col] - feature_value);
                    }
                }

                avg_diff_first /= region_features.n_cols;
                avg_diff_second /= region_features.n_cols;

                float max_allowed_cn_error = 10.0; // Define a threshold for the maximum allowed error
                
                spdlog::debug("Read name: {}", read_name);
                spdlog::debug("Average difference from first cluster: {}", avg_diff_first);
                spdlog::debug("Average difference from second cluster: {}", avg_diff_second);
                spdlog::debug("CN1 diff: {}", cn1_diff);
                spdlog::debug("CN2 diff: {}", cn2_diff);


                if (avg_diff_first <= avg_diff_second && cn1_diff <= max_allowed_cn_error) {
                    spdlog::debug("Adding  read to first cluster");
                    final_clusters[0].emplace_back(read_name, read_intervals);
                } else if (avg_diff_second <= avg_diff_first && cn2_diff <= max_allowed_cn_error) {
                    spdlog::debug("Adding  read to second cluster");
                    final_clusters[1].emplace_back(read_name, read_intervals);
                } else {
                    spdlog::debug("Error too high, adding to noise cluster");
                    noise_cluster.emplace_back(read_name, read_intervals);
                }
            }
        }
    } else if (clusters.size() > 2 && this->manual_reclustering == false && this->somatic == false) {
        if (clusters.size() > 2) {
                // add the other clusters to the noise
            for (size_t i = 2; i< clusters.size();++i){
                auto cluster = clusters[i];
                for (size_t j = 0; j< cluster.size();++j){
                    // print the intervals in the second cluster
                    auto read_name = std::get<0>(reads_intervals[cluster[j]]);
                    auto read_intervals = std::get<1>(reads_intervals[cluster[j]]);
                    noise_cluster.push_back(std::make_pair(read_name, read_intervals));
                }
            }
        }
    }

    // sort the clusters by their size form the largest to the smallest
    std::sort(final_clusters.begin(), final_clusters.end(), [](auto  a, auto b) {
        return a.size() > b.size();
    });

    
    // if there is only one cluster, push and empty vector to the second cluster
    if (final_clusters.size() == 1){
        final_clusters.push_back({});
    }

    
    if (this->somatic == false ){
        final_clusters.resize(2);
    }
    return;
}






