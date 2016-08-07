// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "malformed_file_error.hpp"

#include <utility>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <boost/optional.hpp>
#include <boost/filesystem/operations.hpp>

namespace octopus {

MalformedFileError::MalformedFileError(Path name)
: name_ {std::move(name)}
{}

MalformedFileError::MalformedFileError(Path name, std::string required_type)
: name_ {std::move(name)}, valid_types_ {std::move(required_type)} {}

MalformedFileError::MalformedFileError(Path name, std::vector<std::string> valid_types)
: name_ {std::move(name)}, valid_types_ {std::move(valid_types)} {}

boost::optional<std::string> get_type(const boost::filesystem::path& file)
{
    auto extension = file.extension().string();
    
    if (extension.empty()) return boost::none;
    
    extension.erase(extension.begin());
    
    if (extension == "bam") {
        return std::string {"bam"};
    } else if (extension == "cram") {
        return std::string {"cram"};
    } else if (extension == "bai") {
        return std::string {"bam index"};
    } else if (extension == "crai") {
        return std::string {"cram index"};
    } else if (extension == "fa") {
        return std::string {"fasta"};
    } else if (extension == "fasta") {
        return std::string {"fasta"};
    } else if (extension == "fai") {
        return std::string {"fasta index"};
    } else if (extension == "vcf") {
        return std::string {"vcf"};
    } else if (extension == "bcf") {
        return std::string {"bcf"};
    } else {
        return boost::none;
    }
}

std::string MalformedFileError::do_why() const
{
    std::ostringstream ss {};
    
    const auto type = get_type(name_);
    
    ss << "the ";
    
    if (type) {
        ss << *type << ' ';
    }
    
    ss << "you specified " << name_ << ' ';
    
    if (valid_types_.empty()) {
        ss << "is malformed or currupted";
    } else if (valid_types_.size() == 1) {
        ss << "is not a valid " << valid_types_.front() << " file";
    } else if (valid_types_.size() == 2) {
        ss << "is not a valid " << valid_types_.front() << " or " << valid_types_.back() << " file";
    } else {
        ss << "is not a valid format (from: ";
        std::copy(std::cbegin(valid_types_), std::prev(std::cend(valid_types_)),
                  std::ostream_iterator<std::string> {ss, "; "});
        ss << valid_types_.back() << ')';
    }
    
    return ss.str();
}

std::string MalformedFileError::do_help() const
{
    if (!valid_types_.empty()) {
        return "check the file is not currupted";
    }
    return "check you did not mistake the command line option";
}

} // namespace octopus
