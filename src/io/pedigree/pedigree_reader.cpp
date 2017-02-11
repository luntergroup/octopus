// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree_reader.hpp"

#include <deque>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>

#include <boost/optional.hpp>

#include "exceptions/missing_file_error.hpp"
#include "exceptions/malformed_file_error.hpp"
#include "utils/string_utils.hpp"

namespace octopus { namespace io {

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

PedRecord parse_ped_line(const std::string& line)
{
    static const std::string delims {' ', '\t'};
    auto tokens = utils::split(line, delims);
    return {tokens[0], tokens[1], tokens[4], tokens[2], tokens[3]};
}

bool is_founder(const PedRecord& record)
{
    return record.mother == "0" || record.father == "0";
}

Pedigree::Member::Sex to_pedigree_sex(const std::string& ped_sex)
{
    using Sex = Pedigree::Member::Sex;
    if (ped_sex == "0") {
        return Sex::male;
    } else {
        return Sex::female;
    }
}

} // namespace

Pedigree read_pedigree(const boost::filesystem::path& ped_file)
{
    auto ped_stream = open_ped(ped_file);
    std::deque<PedRecord> records {};
    std::transform(std::istream_iterator<Line> {ped_stream}, std::istream_iterator<Line> {},
                   std::back_inserter(records), [] (const auto& line) { return parse_ped_line(line); });
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
    std::cout << result.size() << std::endl;
    return result;
}
    
} // namespace io
} // namespace octopus
