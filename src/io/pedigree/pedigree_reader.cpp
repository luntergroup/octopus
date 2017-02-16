// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree_reader.hpp"

#include <deque>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <stdexcept>

#include <boost/optional.hpp>
#include <boost/filesystem/operations.hpp>

#include "exceptions/missing_file_error.hpp"
#include "exceptions/malformed_file_error.hpp"
#include "utils/string_utils.hpp"

namespace octopus { namespace io {

class MissingPedigreeFile : public MissingFileError
{
    std::string do_where() const override { return "read_pedigree"; }
public:
    MissingPedigreeFile(boost::filesystem::path p) : MissingFileError {std::move(p), "pedigree"} {};
};

class MalformedPED : public MalformedFileError
{
    std::string do_where() const override { return "read_pedigree"; }
    std::string do_help() const override { return "refer to the latest PED specification"; }
public:
    MalformedPED(boost::filesystem::path file) : MalformedFileError {std::move(file), "ped"} {}
};

namespace {

auto open_ped(const boost::filesystem::path& ped_file)
{
    return std::ifstream {ped_file.string()};
}

struct Line
{
    std::string line_data;
    operator std::string() const { return line_data; }
};

std::istream& operator>>(std::istream& is, Line& data)
{
    std::getline(is, data.line_data);
    if (!data.line_data.empty() && data.line_data.back() == '\r') {
        data.line_data.pop_back();
    }
    return is;
}

struct PedRecord
{
    std::string ped_id, sample, sex, mother, father;
};

bool is_valid_ped_record(const std::vector<std::string>& fields)
{
    return fields.size() >= 5 && (fields[4] == "1" || fields[4] == "2");
}

PedRecord parse_ped_line(const std::string& line)
{
    static const std::string delims {' ', '\t'};
    auto fields = utils::split(line, delims);
    if (!is_valid_ped_record(fields)) {
        throw std::runtime_error {"Malformed PED record"};
    }
    return {fields[0], fields[1], fields[4], fields[2], fields[3]};
}

bool is_founder(const PedRecord& record)
{
    return record.mother == "0" || record.father == "0";
}

Pedigree::Member::Sex to_pedigree_sex(const std::string& ped_sex)
{
    using Sex = Pedigree::Member::Sex;
    if (ped_sex == "1") {
        return Sex::male;
    } else {
        return Sex::female;
    }
}

} // namespace

Pedigree read_pedigree(const boost::filesystem::path& ped_file)
{
    if (!boost::filesystem::exists(ped_file)) {
        throw MissingPedigreeFile {ped_file};
    }
    auto ped_stream = open_ped(ped_file);
    std::deque<PedRecord> records {};
    try {
        std::transform(std::istream_iterator<Line> {ped_stream}, std::istream_iterator<Line> {},
                       std::back_inserter(records), [] (const auto& line) { return parse_ped_line(line); });
    } catch (const std::runtime_error& e) {
        throw MalformedPED {ped_file};
    }
    const auto first_descendant = std::partition(std::begin(records), std::end(records), is_founder);
    Pedigree result {records.size()};
    std::for_each(std::begin(records), first_descendant,
                  [&result] (PedRecord& record) {
                      result.add_founder({record.sample, to_pedigree_sex(record.sex)});
                  });
    std::for_each(first_descendant, std::end(records),
                  [&result] (PedRecord& record) {
                      result.add_descendant({record.sample, to_pedigree_sex(record.sex)}, record.mother, record.father);
                  });
    return result;
}
    
} // namespace io
} // namespace octopus
