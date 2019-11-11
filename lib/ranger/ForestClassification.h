/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef FORESTCLASSIFICATION_H_
#define FORESTCLASSIFICATION_H_

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "globals.h"
#include "Forest.h"

namespace ranger {

class ForestClassification: public Forest {
public:
  ForestClassification() = default;

  ForestClassification(const ForestClassification&) = delete;
  ForestClassification& operator=(const ForestClassification&) = delete;

  virtual ~ForestClassification() override = default;

  void loadForest(size_t num_trees, std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<double>& class_values, std::vector<bool>& is_ordered_variable);

  const std::vector<double>& getClassValues() const {
    return class_values;
  }

  void setClassWeights(std::vector<double>& class_weights) {
    this->class_weights = class_weights;
  }

private:
  // Classes of the dependent variable and classIDs for responses
  std::vector<double> class_values;
  std::vector<uint> response_classIDs;
  std::vector<std::vector<size_t>> sampleIDs_per_class;

  // Splitting weights
  std::vector<double> class_weights;

  // Table with classifications and true classes
  std::map<std::pair<double, double>, size_t> classification_table;
  
  void initInternal() override;
  void growInternal() override;
  void allocatePredictMemory() override;
  void predictInternal(size_t sample_idx) override;
  void computePredictionErrorInternal() override;
  void writeOutputInternal() const override;
  void writeConfusionFile() const override;
  void writePredictionFile() const override;
  void saveToFileInternal(std::ofstream& outfile) const override;
  void loadFromFileInternal(std::ifstream& infile) override;

  double getTreePrediction(size_t tree_idx, size_t sample_idx) const;
  size_t getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const;
};

} // namespace ranger

#endif /* FORESTCLASSIFICATION_H_ */
