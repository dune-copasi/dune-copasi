// Enable latex render
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

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
          to: '/docs/category/docs',
          label: 'Docs',
        },
        {
          label: 'Math Model',
          to: '/math_model',
        },
        {
          to: 'blog',
          label: 'Blog',
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
      links: [
        {
          title: 'Learn',
          items: [
            {
              label: 'About',
              to: '/about',
            },
            {
              label: 'Math Model',
              to: '/math_model',
            },
            {
              label: 'Config Options',
              to: '/docs/param_tree',
            },
            {
              label: 'Docs',
              to: '/docs/category/docs',
            },
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'Community',
              to: '/community',
            },
            {
              label: 'Participate',
              to: '/participate',
            },
          ],
        },
        {
          title: 'More',
          items: [
            {
              label: 'GitLab',
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/dune-copasi/dune-copasi',
            },
            {
              label: 'Docker Registry',
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/container_registry',
            },
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
          title: 'Legal',
          items: [
            {
              href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/blob/master/LICENSE.md',
              label: 'Licence',
            },
            {
              label: 'Impressum',
              href: 'https://www.uni-heidelberg.de/impressum.html',
            },
            {
              label: 'Datenschutzerklärung',
              href: 'https://www.uni-heidelberg.de/datenschutzerklaerung_web.html',
            },
            {
              label: 'Haftungsausschluss',
              href: 'https://www.uni-heidelberg.de/haftungsausschluss_web.html',
            },
          ],
        },
      ],

      logo: {
        href: 'https://www.bmbf.de',
        alt: 'BMBF Logo',
        src: 'img/BMBF_logo.svg',
      },
      copyright: `Research funded by the German federal Ministry of Education and Research (BMBF) FKZ 031L0158</br>Copyright © ${new Date().getFullYear()} Ruprecht-Karls-Universität Heidelberg`,
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
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
          showLastUpdateTime: true,
          showLastUpdateAuthor: true,
        },
        pages: {
          path: 'src/pages',
          routeBasePath: '',
          include: ['**/*.{js,jsx,ts,tsx,md,mdx}'],
          exclude: [
            '**/_*.{js,jsx,ts,tsx,md,mdx}',
            '**/_*/**',
            '**/*.test.{js,jsx,ts,tsx}',
            '**/__tests__/**',
          ],
          mdxPageComponent: '@theme/MDXPage',
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
        },
        blog: {
          path: 'blog',
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
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
