// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef unimplemented_feature_error_hpp
#define unimplemented_feature_error_hpp

#include <string>

#include "program_error.hpp"

namespace octopus {

class UnimplementedFeatureError : public ProgramError
{
public:
    UnimplementedFeatureError() = delete;
    UnimplementedFeatureError(std::string feature, std::string where);
    
    virtual ~UnimplementedFeatureError() override = default;

private:
    virtual std::string do_why() const override;
    virtual std::string do_help() const override;
    virtual std::string do_where() const override;
    
    std::string feature_, where_;
};

} // namespace octopus

#endif
