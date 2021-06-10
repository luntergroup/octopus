const organizationName = 'luntergroup'
const projectName = 'octopus'
const branch = 'master'
const repoUrl = `https://github.com/${organizationName}/${projectName}`

/** @type {import('@docusaurus/types').DocusaurusConfig} */
module.exports = {
  title: 'Octopus',
  tagline: 'Haplotype-based variant calling',
  url: 'https://luntergroup.github.io',
  baseUrl: '/octopus/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  organizationName: organizationName,
  projectName: projectName,
  themeConfig: {
    navbar: {
      hideOnScroll: true,
      title: 'Octopus',
      logo: {
        alt: 'Octopus Logo',
        src: 'img/octopus.svg',
        srcDark: 'img/octopus_keytar.svg',
      },
      items: [
        {
          type: 'doc',
          position: 'left',
          docId: 'introduction',
          label: 'Docs',
        },
        {
          type: 'doc',
          docId: 'cli',
          position: 'left',
          label: 'Command line',
        },
        {
          type: 'doc',
          docId: 'publications',
          position: 'left',
          label: 'Publications',
        },
        {to: '/blog', label: 'Blog', position: 'left'},
        // right
        {
          href: repoUrl,
          position: 'right',
          className: 'header-github-link',
          'aria-label': 'GitHub repository',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Installation',
              to: '/docs/installation',
            },
            {
              label: 'Command Line',
              to: '/docs/cli',
            },
            {
              label: 'Case Studies',
              to: '/docs/tutorials/germline',
            },
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'Biostars',
              href: 'https://www.biostars.org/tag/octopus',
            },
            {
              label: 'Gitter',
              href: 'https://gitter.im/octopus-caller/Lobby',
            },
            {
              label: 'Twitter',
              href: 'https://twitter.com/octopus',
            },
          ],
        },
        {
          title: 'More',
          items: [
            {
              label: 'Blog',
              to: '/blog',
            },
            {
              label: 'GitHub',
              href: repoUrl,
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} Daniel Cooke. Built with Docusaurus.`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/${organizationName}/${projectName}/edit/${branch}/website/',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/${organizationName}/${projectName}/edit/${branch}/website/blog/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
