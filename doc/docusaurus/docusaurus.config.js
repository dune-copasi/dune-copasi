module.exports = {
  title: 'Dune Copasi',
  tagline: 'Solver for reaction-diffusion systems in multiple compartments',
  url: 'https://dune-copasi.github.io',
  baseUrl: '/',
  onBrokenLinks: 'throw',
  favicon: 'img/favicon.ico',
  organizationName: 'dune-copasi',
  projectName: 'dune-copasi',
  themeConfig: {
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
          label: 'Docs',
          position: 'left',
        },
        // {to: 'blog', label: 'Blog', position: 'left'},
        {
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi',
          label: 'GitLab',
          position: 'right',
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
              label: 'Config File',
              to: 'docs/',
            },
          ],
        },
      ],
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          homePageId: 'config_intro',
          sidebarPath: require.resolve('./sidebars.js'),
          editUrl:
          'https://gitlab.dune-project.org/copasi/dune-copasi', // TODO
        },
        blog: {
          showReadingTime: true,
          editUrl:
          'https://gitlab.dune-project.org/copasi/dune-copasi', // TODO
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
