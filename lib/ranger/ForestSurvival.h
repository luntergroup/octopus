/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef FORESTSURVIVAL_H_
#define FORESTSURVIVAL_H_

#include <iostream>
#include <vector>

#include <boost/iostreams/filtering_stream.hpp>

#include "globals.h"
#include "Forest.h"
#include "TreeSurvival.h"

namespace ranger {

class ForestSurvival: public Forest {
public:
  ForestSurvival() = default;

  ForestSurvival(const ForestSurvival&) = delete;
  ForestSurvival& operator=(const ForestSurvival&) = delete;

  virtual ~ForestSurvival() override = default;

  void loadForest(size_t num_trees, std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<std::vector<std::vector<double>> >& forest_chf, std::vector<double>& unique_timepoints,
      std::vector<bool>& is_ordered_variable);

  std::vector<std::vector<std::vector<double>>> getChf() const;

  const std::vector<double>& getUniqueTimepoints() const {
    return unique_timepoints;
  }

private:
  std::vector<double> unique_timepoints;
  std::vector<size_t> response_timepointIDs;
  
  void initInternal() override;
  void growInternal() override;
  void allocatePredictMemory() override;
  void predictInternal(size_t sample_idx) override;
  void computePredictionErrorInternal() override;
  void writeOutputInternal() const override;
  void writeConfusionFile() const override;
  void writePredictionFile() const override;
  void saveToFileInternal(std::ofstream& outfile) const override;
  void loadFromFileInternal(boost::iostreams::filtering_istream& infile) override;

  const std::vector<double>& getTreePrediction(size_t tree_idx, size_t sample_idx) const;
  size_t getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const;
};

} // namespace ranger

#endif /* FORESTSURVIVAL_H_ */
