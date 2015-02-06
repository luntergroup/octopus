//
//  variant_prior_model.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_prior_model__
#define __Octopus__variant_prior_model__

class VariantPriorModel
{
public:
    virtual double get_prior_probability() const noexcept;
};

class SnpPriorModel : public VariantPriorModel
{
    double get_prior_probability() const noexcept override;
};

class MnpPriorModel : public VariantPriorModel
{
    double get_prior_probability() const noexcept override;
};

class InsertionPriorModel : public VariantPriorModel
{
    double get_prior_probability() const noexcept override;
};

class DeletionPriorModel : public VariantPriorModel
{
    double get_prior_probability() const noexcept override;
};

#endif /* defined(__Octopus__variant_prior_model__) */
