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
      href: 'https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.css',
      type: 'text/css',
      integrity: 'sha384-zB1R0rpPzHqg7Kpt0Aljp8JPLqbXI3bhnPWROx27a9N0Ll6ZP/+DiW/UqRcLbRjq',
      crossorigin: 'anonymous',
    },
  ],

  themeConfig: {
    colorMode: {
      defaultMode: 'light',
      disableSwitch: true,
    },
    prism: {
      additionalLanguages: ['ini'],
    },
    announcementBar: {
      id: 'supportus',
      content:
        '⭐️ If you like DuneCopasi, give it a star on <a target="_blank" rel="noopener noreferrer" href="https://gitlab.dune-project.org/copasi/dune-copasi">GitLab</a> or <a target="_blank" rel="noopener noreferrer" href="https://github.com/dune-copasi/dune-copasi">GitHub</a>! ⭐️',
    },
    navbar: {
      title: 'DuneCopasi',
      logo: {
        alt: 'dune-copasi',
        src: 'img/logo.svg',
      },
      items: [
        {
          to: 'docs/',
          activeBasePath: 'docs',
          label: 'Learn',
          position: 'left',
        },
        {
          position: 'left',
          label: 'Help',
          items: [
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/issues',
              label: 'Bug tracker',
              position: 'left',
            },
          ],
        },
        // {to: 'blog', label: 'Blog', position: 'left'},
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
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          homePageId: 'about',
          sidebarPath: require.resolve('./sidebars.js'),
          editUrl:
          'https://gitlab.dune-project.org/copasi/dune-copasi/-/edit/master/doc/docusaurus/',
          remarkPlugins: [math],
          rehypePlugins: [katex],
          showLastUpdateTime: true,
          showLastUpdateAuthor: true,
        },
        blog: {
          showReadingTime: true,
          editUrl:
          'https://gitlab.dune-project.org/copasi/dune-copasi/-/edit/master/doc/docusaurus/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
