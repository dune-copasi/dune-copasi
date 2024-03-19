module.exports = {
  docs: [
    {
      type: 'category',
      label: 'Information',
      link: {
        type: 'generated-index',
        slug: '/category/info',
      },
      items: [
        'about',
        'math_model',
        'community',
        'participate',
        {
          type: 'link',
          label: 'Releases',
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/releases',
          description: "List of stable releases",
        },
        {
          type: 'link',
          label: 'License',
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/blob/master/LICENSE.md',
          description: "BSD 2-Clause \"Simplified\" License",
        },
        {
          type: 'link',
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/blob/master/CHANGELOG.md',
          label: 'Changelog',
          description: "Notable changes on the project",
        },
        {
          type: 'link',
          href: 'https://gitlab.dune-project.org/copasi/dune-copasi/-/issues?label_name%5B%5D=Bug',
          label: 'Bug tracker',
          description: "Known issues",
        },
      ],
    },
    {
      type: 'category',
      label: 'Documentation',
      link: {
        type: 'generated-index',
        slug: '/category/docs',
      },
      items: [
        'install_use',
        {
          type: 'category',
          label: 'Configuration Options',
          link: {
            type: 'generated-index',
            slug: '/category/docs/config',
          },
          items: [
            'ini_file',
            'param_tree',
            'input_data',
            'math_expr',
          ]
        },
        'output',
        'api',
      ]
    }
  ],
};