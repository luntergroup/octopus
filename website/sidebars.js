/**
 * Creating a sidebar enables you to:
 - create an ordered group of docs
 - render a sidebar for each doc of that group
 - provide next/previous navigation

 The sidebars can be generated from the filesystem, or explicitly defined here.

 Create as many sidebars as you want.
 */

module.exports = {
  docs: [
    {
      type: 'category',
      label: 'Getting Started',
      collapsed: false,
      items: [
        'introduction',
        'installation',
      ],
    },
    {
      type: 'category',
      label: 'Guides',
      collapsed: true,
      items: [
        'guides/overview',
        'guides/preprocessing',
        'guides/discovery',
        'guides/haplotypes',
        {
          type: 'category',
          label: 'Calling Models',
          collapsed: true,
          items: [
            'guides/models/individual',
            'guides/models/trio',
            'guides/models/population',
            'guides/models/cancer',
            'guides/models/polyclone',
            'guides/models/cell',
          ]
        },
        'guides/errorModels',
        'guides/phasing',
        {
          Filtering: [
            'guides/filtering/annotations',
            'guides/filtering/forest',
            'guides/filtering/thresholds',
          ]
        },
        'guides/bamout',
        {
          Advanced: [
            'guides/advanced/targeted',
            'guides/advanced/threading',
            'guides/advanced/memory',
            'guides/advanced/vcf',
            'guides/advanced/configs',
          ]
        },
      ],
    },
    {
      type: 'category',
      label: 'Case studies',
      collapsed: true,
      items: [
        'tutorials/germline',
        'tutorials/mitochondria'
      ],
    }
  ],
};
