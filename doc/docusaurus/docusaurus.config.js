// Enable latex render
const math = require('remark-math')
const katex = require('rehype-katex')

module.exports = {
  title: 'Dune Copasi',
  tagline: 'Solver for reaction-diffusion systems in multiple compartments',
  url: 'https://dune-copasi.netlify.app/',
  baseUrl: '/',
  onBrokenLinks: 'throw',
  favicon: 'img/favicon.ico',
  organizationName: 'dune-copasi',
  projectName: 'dune-copasi',

  // CSS for latex
  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css',
      type: 'text/css',
      integrity: 'sha384-AfEj0r4/OFrOo5t7NnNe46zW/tFgW6x/bCJG8FqQCEo3+Aro6EYUG4+cU+KJWu/X',
      crossorigin: 'anonymous',
    },
  ],

  themeConfig: {
    colorMode: {
      defaultMode: 'light',
      disableSwitch: true,
    },
    prism: {
      additionalLanguages: ['ini','cmake'],
    },
    navbar: {
      title: 'DuneCopasi',
      logo: {
        alt: 'dune-copasi',
        src: 'img/logo.svg',
      },
      items: [
        {
          to: 'docs/math_model',
          label: 'Learn',
          position: 'left',
        },
        {
          type: 'docsVersionDropdown',
          position: 'left',
          dropdownActiveClassDisabled: true,
        },
        {
          position: 'right',
          label: 'Help',
          items: [
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/releases',
              label: 'Releases',
            },
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/blob/master/CHANGELOG.md',
              label: 'Changelog',
            },
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/issues?label_name%5B%5D=Bug',
              label: 'Bug tracker',
            },
          ],
        },
        {
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi',
          className: 'header-gitlab-link',
          position: 'right',
          'aria-label': 'GitLab repository',
        },
        {
          href: 'https://github.com/dune-copasi/dune-copasi',
          position: 'right',
          className: 'header-github-link',
          'aria-label': 'GitHub repository',
        },
        {
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi/container_registry',
          position: 'right',
          className: 'header-docker-link',
          'aria-label': 'Docker registry',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [],
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          path: 'docs',
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/edit/master/doc/docusaurus/',
          remarkPlugins: [ math ],
          rehypePlugins: [ katex ],
          showLastUpdateTime: true,
          showLastUpdateAuthor: true,
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
