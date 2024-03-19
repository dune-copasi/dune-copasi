module.exports = {
  docs: [
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