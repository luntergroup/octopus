// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mock_variant_generator_hpp
#define mock_variant_generator_hpp

#include <core/tools/variantgenerator/variant_generator.hpp>

namespace octopus { namespace test {

    class MockVariantGenerator : public VariantGenerator
    {
    public:
        MockVariantGenerator() = delete;
        
        MockVariantGenerator(const MockVariantGenerator&)            = default;
        MockVariantGenerator& operator=(const MockVariantGenerator&) = default;
        MockVariantGenerator(MockVariantGenerator&&)                 = default;
        MockVariantGenerator& operator=(MockVariantGenerator&&)      = default;
        
        ~MockVariantGenerator() override = default;
        
    private:
        
    };

} // namespace test
} // namespace octopus

#endif
