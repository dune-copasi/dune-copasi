// Enable latex render
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

const copyright = `Copyright © ${new Date().getFullYear()} Ruprecht-Karls-Universität Heidelberg`
const impressum = `<a href="https://www.uni-heidelberg.de/impressum.html" rel="noreferrer">Impressum</a>`
const datenschutz = `<a href="https://www.uni-heidelberg.de/datenschutzerklaerung_web.html" rel="noreferrer">Datenschutzerklärung</a>`
const disclaimer = `<a href="https://www.uni-heidelberg.de/haftungsausschluss_web.html" rel="noreferrer">Haftungsausschluss</a>`


module.exports = {
  title: 'Dune Copasi',
  tagline: 'Solver for reaction-diffusion systems in multiple compartments',
  url: 'https://dune-copasi.netlify.app/',
  baseUrl: '/',
  onBrokenLinks: 'throw',
  favicon: 'img/favicon.ico',
  projectName: 'dune-copasi',

  // CSS for latex
  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.13.24/dist/katex.min.css',
      type: 'text/css',
      integrity: 'sha384-odtC+0UGzzFL/6PNoE8rX/SPcQDXBJ+uRepguP4QkPCm2LBxH3FA3y+fKSiJ+AmM',
      crossorigin: 'anonymous',
    },
  ],

  themeConfig: {
    colorMode: {
      defaultMode: 'light',
      disableSwitch: true,
    },
    prism: {
      additionalLanguages: ['ini', 'cmake', 'bash'],
      magicComments: [
        // First one: default highlighted line by docusaurus
        {
          className: 'theme-code-block-highlighted-line',
          line: '%highlight',
          block: { start: '%highlight-start', end: '%highlight-end' },
        },
        // Custom highlighted line
        {
          className: 'code-block-highlighted-error-line',
          line: '%highlight-error',
          block: { start: '%highlight-start-error', end: '%highlight-end-error' },
        },
      ],
    },
    navbar: {
      title: 'DuneCopasi',
      logo: {
        alt: 'dune-copasi',
        src: 'img/logo.svg',
      },
      items: [
        {
          to: 'docs/category/info',
          label: 'Info',
          position: 'left',
        },
        {
          to: 'docs/category/docs',
          label: 'Docs',
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
              label: '🚀 Releases',
            },
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/blob/master/LICENSE.md',
              label: '⚖️ Licence',
            },
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/blob/master/CHANGELOG.md',
              label: '🚧 Changelog',
            },
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/issues?label_name%5B%5D=Bug',
              label: '🐛 Bug tracker',
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
      style: 'light',
      links: [],
      logo: {
        alt: 'Meta Open Source Logo',
        src: 'img/BMBF_logo.svg',
      },
      // copyright: `Research funded by the German federal Ministry of Education and Research (BMBF) FKZ 031L0158</br>${copyright} | ${impressum} | ${datenschutz} | ${disclaimer}`,
      copyright: `Research funded by the German federal Ministry of Education and Research (BMBF) FKZ 031L0158`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          path: 'docs',
          sidebarPath: require.resolve('./sidebars.js'),
          editUrl: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/edit/master/doc/docusaurus/',
          remarkPlugins: [ remarkMath ],
          rehypePlugins: [ rehypeKatex ],
          showLastUpdateTime: true,
          showLastUpdateAuthor: true,
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
  plugins: [
    require.resolve("@cmfcmf/docusaurus-search-local"),
    process.env.NODE_ENV === 'production' && '@docusaurus/plugin-debug',
  ].filter(Boolean),
};
