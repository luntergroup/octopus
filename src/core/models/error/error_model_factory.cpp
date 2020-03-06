// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_model_factory.hpp"

#include <array>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>

#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>

#include "utils/string_utils.hpp"
#include "exceptions/user_error.hpp"
#include "exceptions/malformed_file_error.hpp"

#include "basic_repeat_based_indel_error_model.hpp"
#include "repeat_based_snv_error_model.hpp"
#include "custom_repeat_based_indel_error_model.hpp"

namespace octopus {

static constexpr std::array<LibraryPreparation, 3> libraries {
    LibraryPreparation::pcr, LibraryPreparation::pcr_free, LibraryPreparation::tenx
};
static constexpr std::array<Sequencer, 8> sequencers {
    Sequencer::hiseq_2000, Sequencer::hiseq_2500, Sequencer::hiseq_4000,
    Sequencer::xten, Sequencer::novaseq, Sequencer::bgiseq_500,
    Sequencer::pacbio, Sequencer::pacbio_css
};

std::ostream& operator<<(std::ostream& out, const LibraryPreparation& library)
{
    switch (library) {
        case LibraryPreparation::pcr:
            out << "PCR";
            break;
        case LibraryPreparation::pcr_free:
            out << "PCR-free";
            break;
        case LibraryPreparation::tenx:
            out << "10X";
            break;
    }
    return out;
}

template <typename Range>
std::ostream& join(const Range& range, std::ostream& os, const std::string& delim = " ")
{
    if (std::cbegin(range) != std::cend(range)) {
        using T = typename std::iterator_traits<decltype(std::cbegin(range))>::value_type;
        std::copy(std::cbegin(range), std::prev(std::cend(range)), std::ostream_iterator<T> {os, delim.c_str()});
        os << *std::prev(std::cend(range));
    }
    return os;
}

class UnknownLibraryPreparation : public UserError
{
    std::string name_;
    
    std::string do_where() const override
    {
        return "operator>>(std::istream&, LibraryPreparation&)";
    }
    std::string do_why() const override
    {
        return "The library preparation name " + name_ + " is unknown";
    }
    std::string do_help() const override
    {
        std::ostringstream ss {};
        ss << "Choose a valid library preparation name [";
        join(libraries, ss, ", ");
        ss << "]";
        return ss.str();
    }
public:
    UnknownLibraryPreparation(std::string name) : name_ {std::move(name)} {}
};

std::istream& operator>>(std::istream& in, LibraryPreparation& result)
{
    std::string token;
    in >> token;
    utils::capitalise(token);
    if (token == "PCR")
        result = LibraryPreparation::pcr;
    else if (token == "PCR-FREE" || token == "PCRF")
        result = LibraryPreparation::pcr_free;
    else if (token == "10X")
        result = LibraryPreparation::tenx;
    else throw UnknownLibraryPreparation {token};
    return in;
}

std::ostream& operator<<(std::ostream& out, const Sequencer& sequencer)
{
    switch (sequencer) {
        case Sequencer::hiseq_2000:
            out << "HiSeq-2000";
            break;
        case Sequencer::hiseq_2500:
            out << "HiSeq-2500";
            break;
        case Sequencer::hiseq_4000:
            out << "HiSeq-4000";
            break;
        case Sequencer::xten:
            out << "X10";
            break;
        case Sequencer::novaseq:
            out << "NovaSeq";
            break;
        case Sequencer::bgiseq_500:
            out << "BGISEQ-500";
            break;
        case Sequencer::pacbio:
            out << "PacBio";
            break;
        case Sequencer::pacbio_css:
            out << "PacBioCSS";
            break;
    }
    return out;
}

class UnknownSequencer : public UserError
{
    std::string name_;
    
    std::string do_where() const override
    {
        return "operator>>(std::istream&, Sequencer&)";
    }
    std::string do_why() const override
    {
        return "The sequencer name " + name_ + " is unknown";
    }
    std::string do_help() const override
    {
        std::ostringstream ss {};
        ss << "Choose a valid sequencer name [";
        join(sequencers, ss, ", ");
        ss << "]";
        return ss.str();
    }
public:
    UnknownSequencer(std::string name) : name_ {std::move(name)} {}
};

std::istream& operator>>(std::istream& in, Sequencer& result)
{
    std::string token;
    in >> token;
    utils::capitalise(token);
    if (token == "HISEQ-2000")
        result = Sequencer::hiseq_2000;
    else if (token == "HISEQ-2500")
        result = Sequencer::hiseq_2500;
    else if (token == "HISEQ-4000")
        result = Sequencer::hiseq_4000;
    else if (token == "X10")
        result = Sequencer::xten;
    else if (token == "NOVASEQ")
        result = Sequencer::novaseq;
    else if (token == "BGISEQ-500")
        result = Sequencer::bgiseq_500;
    else if (token == "PACBIO")
        result = Sequencer::pacbio;
    else if (token == "PACBIOCSS")
        result = Sequencer::pacbio_css;
    else throw UnknownSequencer {token};
    return in;
}

LibraryPreparation to_library(const std::string& name)
{
    LibraryPreparation result;
    std::istringstream ss {name};
    ss >> result;
    return result;
}

Sequencer to_sequencer(const std::string& name)
{
    Sequencer result;
    std::istringstream ss {name};
    ss >> result;
    return result;
}

struct ModelConfigHash
{
    std::size_t operator()(const ModelConfig& config) const
    {
        using boost::hash_combine;
        std::size_t seed {};
        hash_combine(seed, config.library);
        hash_combine(seed, config.sequencer);
        return seed;
    }
};

bool operator==(const ModelConfig& lhs, const ModelConfig& rhs) noexcept
{
    return lhs.library == rhs.library && lhs.sequencer == rhs.sequencer;
}

using RepeatBasedIndelModelParameterMap = std::unordered_map<ModelConfig, BasicRepeatBasedIndelErrorModel::Parameters, ModelConfigHash>;

static const RepeatBasedIndelModelParameterMap builtin_indel_models {{
    {
        {LibraryPreparation::pcr_free, Sequencer::hiseq_2000},
        {
            {45,45,43,43,41,38,35,32,29,25,21,20,19,18,17,17,16,16,15,14,14,13,12,12,11,10,9,9,8,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5},
            {45,45,45,41,39,34,30,24,21,18,15,13,12,10,8,7,7,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3},
            {45,45,42,40,35,29,26,24,22,21,20,19,18,18,17,17,16,16,15,15,15,14,13,13,12,12,11,11,10,10,9,9,9,7,7,7,6,6,5,4,4,4,4,4,4,4,4,4,3},
            {45,45,40,36,30,28,26,25,23,22,22,22,21,21,20,20,20,18,17,16,14,14,14,14,12,11,11,11,10,10,10,7,7,7,4,4,4,4,4,4,4,3}
        }
    },
    {
        {LibraryPreparation::pcr_free, Sequencer::hiseq_2500},
        {
            {45,45,43,43,41,38,35,32,29,25,21,20,19,18,17,17,16,16,15,14,14,13,12,12,11,10,9,9,8,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5},
            {45,45,45,41,39,34,30,24,21,18,15,13,12,10,8,7,7,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3},
            {45,45,42,40,35,29,26,24,22,21,20,19,18,18,17,17,16,16,15,15,15,14,13,13,12,12,11,11,10,10,9,9,9,7,7,7,6,6,5,4,4,4,4,4,4,4,4,4,3},
            {45,45,40,36,30,28,26,25,23,22,22,22,21,21,20,20,20,18,17,16,14,14,14,14,12,11,11,11,10,10,10,7,7,7,4,4,4,4,4,4,4,3}
        }
    },
    {
     	{LibraryPreparation::pcr_free, Sequencer::hiseq_4000},
        {
            {45,45,43,43,41,38,35,32,29,25,21,20,19,18,17,17,16,16,15,14,14,13,12,12,11,10,9,9,8,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5},
            {45,45,45,41,39,34,30,24,21,18,15,13,12,10,8,7,7,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3},
            {45,45,42,40,35,29,26,24,22,21,20,19,18,18,17,17,16,16,15,15,15,14,13,13,12,12,11,11,10,10,9,9,9,7,7,7,6,6,5,4,4,4,4,4,4,4,4,4,3},
            {45,45,40,36,30,28,26,25,23,22,22,22,21,21,20,20,20,18,17,16,14,14,14,14,12,11,11,11,10,10,10,7,7,7,4,4,4,4,4,4,4,3}
        }
    },
    {
     	{LibraryPreparation::pcr_free, Sequencer::xten},
        {
            {45,45,47,44,43,39,37,33,29,26,24,22,21,21,20,19,18,18,17,16,15,14,14,13,13,13,12,12,11,11,10,10,10,9,9,9,9,9,9,9,9,9,9,9,8},
            {45,45,47,45,42,38,33,28,24,21,19,18,16,15,13,12,11,11,11,10,10,9,9,9,9,9,9,9,8,8,8,8,8,8,8,6},
            {45,45,43,40,35,31,28,26,24,23,22,21,20,20,19,19,18,17,17,17,16,15,15,14,13,13,12,11,11,11,10,10,10,8,8,8,6,6,5},
            {45,45,42,36,32,29,27,26,24,23,23,22,21,21,20,19,18,17,15,14,13,13,13,11,11,11,11,10,7,7,7,5,5,5,5,5,4}
        }
    },
    {
     	{LibraryPreparation::pcr_free, Sequencer::novaseq},
        {
            {45,45,43,42,41,38,35,32,29,26,24,22,21,21,20,19,18,17,16,16,15,14,13,12,11,11,10,9,8,7,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3},
            {45,45,43,40,37,33,28,21,19,17,15,13,12,11,8,7,7,6,6,6,5,5,5,5,3},
            {45,45,39,38,33,29,26,25,23,22,21,20,20,19,18,18,17,17,16,15,15,14,14,13,12,11,10,9,9,9,8,8,8,8,8,6,5,5,5,4,4,4,4,4,4,3},
            {45,45,39,34,30,27,25,25,22,22,21,21,21,21,19,19,19,16,13,13,13,11,11,11,11,11,10,9,8,7,6,4,4,4,4,4,4,4,3}
        }
    },
    {
     	{LibraryPreparation::pcr_free, Sequencer::bgiseq_500},
        {
            {45,45,47,44,43,39,37,33,29,26,24,22,21,21,20,19,18,18,17,16,15,14,14,13,13,13,12,12,11,11,10,10,10,9,9,9,9,9,9,9,9,9,9,9,8},
            {45,45,47,45,42,38,33,28,24,21,19,18,16,15,13,12,11,11,11,10,10,9,9,9,9,9,9,9,8,8,8,8,8,8,8,6},
            {45,45,43,40,35,31,28,26,24,23,22,21,20,20,19,19,18,17,17,17,16,15,15,14,13,13,12,11,11,11,10,10,10,8,8,8,6,6,5},
            {45,45,42,36,32,29,27,26,24,23,23,22,21,21,20,19,18,17,15,14,13,13,13,11,11,11,11,10,7,7,7,5,5,5,5,5,4}
        }
    },
    {
        {LibraryPreparation::pcr_free, Sequencer::pacbio},
        {
            {13,13,11,10,9,8,7,7,7,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4},
            {13,13,10,8,7,7,7,7,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4},
            {13,13,8,7,6,6,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3},
            {13,13,7,6,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}
        }
    },
    {
        {LibraryPreparation::pcr_free, Sequencer::pacbio_css},
        {
            {30,30,26,23,20,18,16,14,13,11,11,10,9,9,9,8,8,8,8,7,7,7,7,7,7,7,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,3},
            {30,30,25,21,18,16,14,12,10,9,8,8,6,6,6,6,6,4,4,4,4,3},
            {30,30,23,21,19,16,14,13,11,10,10,9,9,9,9,8,8,7,7,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,4,4,4,4},
            {30,30,21,19,16,14,13,12,11,11,11,11,11,11,9,9,9,7,6,6,6,6,6,6,6,5}
        }
    },
    {
        {LibraryPreparation::pcr, Sequencer::hiseq_2000},
        {
            {45,45,43,41,40,36,34,30,24,20,16,13,12,11,10,10,9,9,8,8,7,7,7,6,6,6,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {45,45,42,40,37,33,27,21,17,15,12,10,9,7,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,38,37,32,26,21,18,16,14,14,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,7,7,7,7,6,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,37,32,26,22,20,19,18,17,17,16,15,15,14,13,13,12,12,12,12,10,10,10,9,9,7,7,7,7,6,6,6,6,4,3}
        }
    },
    {
        {LibraryPreparation::pcr, Sequencer::hiseq_2500},
        {
            {45,45,43,41,40,36,34,30,24,20,16,13,12,11,10,10,9,9,8,8,7,7,7,6,6,6,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {45,45,42,40,37,33,27,21,17,15,12,10,9,7,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,38,37,32,26,21,18,16,14,14,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,7,7,7,7,6,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,37,32,26,22,20,19,18,17,17,16,15,15,14,13,13,12,12,12,12,10,10,10,9,9,7,7,7,7,6,6,6,6,4,3}
        }
    },
    {
     	{LibraryPreparation::pcr, Sequencer::hiseq_4000},
        {
            {45,45,43,41,40,36,34,30,24,20,16,13,12,11,10,10,9,9,8,8,7,7,7,6,6,6,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {45,45,42,40,37,33,27,21,17,15,12,10,9,7,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,38,37,32,26,21,18,16,14,14,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,7,7,7,7,6,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,37,32,26,22,20,19,18,17,17,16,15,15,14,13,13,12,12,12,12,10,10,10,9,9,7,7,7,7,6,6,6,6,4,3}
        }
    },
    {
     	{LibraryPreparation::pcr, Sequencer::xten},
        {
            {60,60,44,42,40,36,34,30,24,20,16,13,12,11,10,9,9,8,8,8,7,7,7,6,6,6,5,5,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3},
            {60,60,42,40,37,35,28,22,18,15,12,10,9,7,6,4,4,5,5,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3},
            {60,60,38,37,33,27,21,18,16,15,14,13,13,12,12,12,11,11,10,10,9,9,9,9,8,7,7,6,6,6,5,5,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3},
            {60,60,38,33,26,22,20,19,18,18,17,16,17,15,14,14,14,14,13,12,12,11,10,10,8,7,7,6,6,5,5,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}
        }
    },
    {
     	{LibraryPreparation::pcr, Sequencer::novaseq},
        {
            {45,45,43,41,40,36,34,30,24,20,16,13,12,11,10,10,9,9,8,8,7,7,7,6,6,6,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {45,45,42,40,37,33,27,21,17,15,12,10,9,7,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,38,37,32,26,21,18,16,14,14,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,7,7,7,7,6,6,6,5,5,4,4,4,4,4,4,4,3},
            {45,45,37,32,26,22,20,19,18,17,17,16,15,15,14,13,13,12,12,12,12,10,10,10,9,9,7,7,7,7,6,6,6,6,4,3}
        }
    },
    {
     	{LibraryPreparation::pcr, Sequencer::bgiseq_500},
        {
            {60,60,49,47,43,39,35,31,25,21,17,14,13,12,11,11,10,10,9,9,9,8,8,8,8,8,7,7,7,7,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,4},
            {60,60,48,45,42,38,32,26,22,17,14,12,10,9,8,7,6,6,6,5,5,5,5,5,4,4,4,4,4,4,4,4,3},
            {60,60,44,42,36,29,22,19,17,15,15,14,14,13,13,13,12,12,12,11, 11,10,10,10,9,9,9,8,8,8,8,8,8,7,7,6,6,6,5,4,4,4,4,3},
            {60,60,41,36,28,23,21,20,19,18,18,17,17,16,15,15,14,13,12,12,12,12,10,9,9,9,9,8,8,7,7,7,7,6,6,6,6,5,5,5,5,4,4,4,4,4,4,4,3}
        }
    },
    {
        {LibraryPreparation::pcr, Sequencer::pacbio},
        {
            {13,13,11,10,9,8,7,7,7,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4},
            {13,13,10,8,7,7,7,7,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4},
            {13,13,8,7,6,6,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3},
            {13,13,7,6,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}
        }
    },
    {
        {LibraryPreparation::pcr, Sequencer::pacbio_css},
        {
            {40,40,31,29,28,24,21,19,17,15,13,12,11,10,10,9,9,8,8,8,7,7,6,6,6,6,5,5,5,5,4},
            {40,40,33,31,28,22,17,13,12,10,9,8,7,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {40,40,30,27,22,18,16,15,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,8,7,7,6,6,6,6,5,5,5,4},
            {40,40,28,25,19,16,15,14,12,12,12,12,11,11,10,10,10,9,9,9,8,7,7,7,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,3}
        }
    },
    {
        {LibraryPreparation::tenx, Sequencer::hiseq_2000},
        {
            {45,45,36,34,30,27,26,24,20,16,13,12,11,10,9,8,8,7,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4},
            {45,45,34,31,28,24,21,18,16,14,12,10,9,8,8,8,7,7,7,7,7,7,7,7,7,7,6,3},
            {45,45,34,33,29,23,18,15,14,13,12,12,11,11,10,10,10,9,9,9,9,8,8,8,8,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5},
            {45,45,32,29,23,19,17,16,15,14,14,13,12,12,11,11,11,11,11,9,9,9,7,7,7,7,6}
        }
    },
    {
        {LibraryPreparation::tenx, Sequencer::hiseq_2500},
        {
            {45,45,37,35,30,28,26,25,21,17,14,12,11,11,10,10,9,9,8,8,8,7,7,7,7,7,7,7,6},
            {45,45,36,33,29,26,22,20,17,15,13,11,10,9,9,9,9,8},
            {45,45,33,32,28,23,18,15,13,12,12,11,11,10,10,9,9,9,8,8,8,8,7,7,7,6,6,6,6,6,5},
            {45,45,31,28,23,19,16,15,14,13,12,12,11,10,10,10,10,10,10,9,9,7,7,6,6,6,5}
        }
    },
    {
     	{LibraryPreparation::tenx, Sequencer::hiseq_4000},
        {
            {45,45,37,35,30,28,26,25,21,17,14,12,11,11,10,10,9,9,8,8,8,7,7,7,7,7,7,7,6},
            {45,45,36,33,29,26,22,20,17,15,13,11,10,9,9,9,9,8},
            {45,45,33,32,28,23,18,15,13,12,12,11,11,10,10,9,9,9,8,8,8,8,7,7,7,6,6,6,6,6,5},
            {45,45,31,28,23,19,16,15,14,13,12,12,11,10,10,10,10,10,10,9,9,7,7,6,6,6,5}
        }
    },
    {
     	{LibraryPreparation::tenx, Sequencer::xten},
        {
            {45,45,31,29,28,24,21,19,17,15,13,12,11,10,10,9,9,8,8,8,7,7,6,6,6,6,5,5,5,5,4},
            {45,45,33,31,28,22,17,13,12,10,9,8,7,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {45,45,30,27,22,18,16,15,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,8,7,7,6,6,6,6,5,5,5,4},
            {45,45,28,25,19,16,15,14,12,12,12,12,11,11,10,10,10,9,9,9,8,7,7,7,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,3}
        }
    },
    {
     	{LibraryPreparation::tenx, Sequencer::novaseq},
        {
            {45,45,31,29,28,24,21,19,17,15,13,12,11,10,10,9,9,8,8,8,7,7,6,6,6,6,5,5,5,5,4},
            {45,45,33,31,28,22,17,13,12,10,9,8,7,6,5,5,5,4,4,4,4,4,4,4,4,4,4,3},
            {45,45,30,27,22,18,16,15,13,13,12,12,11,11,11,10,10,10,9,9,9,8,8,8,7,7,6,6,6,6,5,5,5,4},
            {45,45,28,25,19,16,15,14,12,12,12,12,11,11,10,10,10,9,9,9,8,7,7,7,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,3}
        }
    },
    {
     	{LibraryPreparation::tenx, Sequencer::bgiseq_500},
        {
            {45,45,37,35,30,28,26,25,21,17,14,12,11,11,10,10,9,9,8,8,8,7,7,7,7,7,7,7,6},
            {45,45,36,33,29,26,22,20,17,15,13,11,10,9,9,9,9,8},
            {45,45,34,33,29,23,18,15,14,13,12,12,11,11,10,10,10,9,9,9,9,8,8,8,8,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5},
            {45,45,32,29,23,19,17,16,15,14,14,13,12,12,11,11,11,11,11,9,9,9,7,7,7,7,6}
        }
    }
 }};

BasicRepeatBasedIndelErrorModel::Parameters lookup_builtin_indel_model(const ModelConfig config)
{
    return builtin_indel_models.at(config);
}

bool use_snv_error_model(const ModelConfig config)
{
    return !(config.sequencer != Sequencer::pacbio || config.sequencer != Sequencer::pacbio_css);
}

using RepeatBasedSnvModelParameterMap = std::unordered_map<LibraryPreparation, BasicRepeatBasedSNVErrorModel::Parameters>;

static const RepeatBasedSnvModelParameterMap builtin_snv_models {{
    {
        LibraryPreparation::pcr_free,
     {
        {125,125,60,55,50,30,20,15,12,12,10,10,10,10,8,7,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1},
        {125,125,60,60,52,52,38,38,22,22,17,17,15,15,13,13,10,10,10,10,8,8,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1},
        {125,125,125,55,55,55,40,40,40,25,25,25,19,19,19,11,11,11,9,9,9,7,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1}
     }
    },
    {
        LibraryPreparation::pcr,
    {
        {125,125,60,55,38,23,16,14,11,10,9,8,7,7,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1},
        {125,125,60,60,52,52,38,38,22,22,17,17,15,15,13,13,10,10,10,10,8,8,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1},
        {125,125,125,55,55,55,40,40,40,25,25,25,19,19,19,11,11,11,9,9,9,7,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1}
    }},
    {
      LibraryPreparation::tenx,
    {
      {125,125,60,55,38,23,16,14,11,10,9,8,7,7,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1},
      {125,125,60,60,52,52,38,38,22,22,17,17,15,15,13,13,10,10,10,10,8,8,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1},
      {125,125,125,55,55,55,40,40,40,25,25,25,19,19,19,11,11,11,9,9,9,7,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1}
    }}
}};

boost::optional<BasicRepeatBasedSNVErrorModel::Parameters> lookup_builtin_snv_model(const ModelConfig config)
{
    if (use_snv_error_model(config)) {
        return builtin_snv_models.at(config.library);
    } else {
        return boost::none;
    }
}

class MalformedErrorModelFile : public MalformedFileError
{
    std::string do_where() const override { return "make_indel_error_model"; }
    std::string do_help() const override { return "refer to documentation on custom error models or use provided Python script"; }
public:
    MalformedErrorModelFile(boost::filesystem::path file) : MalformedFileError {std::move(file), "model"} {}
};

std::unique_ptr<SnvErrorModel> make_snv_error_model(const ModelConfig config)
{
    auto model = lookup_builtin_snv_model(config);
    if (model) {
        return std::make_unique<BasicRepeatBasedSNVErrorModel>(std::move(*model));
    } else {
        return nullptr;
    }
}

std::unique_ptr<IndelErrorModel> make_indel_error_model(const ModelConfig config)
{
    return std::make_unique<BasicRepeatBasedIndelErrorModel>(lookup_builtin_indel_model(config));
}

ModelConfig parse_model_config(const std::string& label)
{
    auto result = default_model_config;
    const auto library_end_pos = label.find('.');
    const auto library_name = label.substr(0, library_end_pos);
    if (!library_name.empty()) {
        result.library = to_library(library_name);
    }
    if (library_end_pos != std::string::npos) {
        const auto sequencer_name = label.substr(library_end_pos + 1);
        if (!sequencer_name.empty()) {
            result.sequencer = to_sequencer(sequencer_name);
        }
    }
    return result;
}

ErrorModel make_error_model(const std::string& label)
{
    const auto config = parse_model_config(label);
    return {make_indel_error_model(config), make_snv_error_model(config)};
}

ErrorModel make_error_model(const boost::filesystem::path& model_filename)
{
    std::ifstream model_file {model_filename.string()};
    std::string model_str {static_cast<std::stringstream const&>(std::stringstream() << model_file.rdbuf()).str()};
    auto params = make_penalty_map(std::move(model_str));
    if (!params.open) {
        throw MalformedErrorModelFile {model_filename};
    }
    ErrorModel result {};
    if (params.extend) {
        result.indel = std::make_unique<CustomRepeatBasedIndelErrorModel>(std::move(*params.open), std::move(*params.extend));
    } else {
        result.indel = std::make_unique<CustomRepeatBasedIndelErrorModel>(std::move(*params.open));
    }
    result.snv = make_snv_error_model(default_model_config);
    return result;
}

} // namespace octopus
