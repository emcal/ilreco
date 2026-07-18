import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'ilreco',
  description: 'Island-clustering reconstruction for square-cell calorimeters',
  base: '/ilreco/',
  themeConfig: {
    nav: [
      { text: 'Install', link: '/install' },
      { text: 'Getting started', link: '/getting-started' },
      { text: 'Python', link: '/python-api' },
      { text: 'C++', link: '/cpp-api' },
    ],
    sidebar: [
      {
        text: 'Guide',
        items: [
          { text: 'Overview', link: '/' },
          { text: 'Install', link: '/install' },
          { text: 'Getting started', link: '/getting-started' },
        ],
      },
      {
        text: 'Using the library',
        items: [
          { text: 'C++ API', link: '/cpp-api' },
          { text: 'Python API', link: '/python-api' },
          { text: 'Threading model', link: '/threading' },
          { text: 'Memory model', link: '/memory-model' },
          { text: 'Advanced multithreading', link: '/advanced-mt' },
        ],
      },
      {
        text: 'Reference',
        items: [
          { text: 'Shower profiles', link: '/profiles' },
          { text: 'Development', link: '/development' },
        ],
      },
    ],
    socialLinks: [
      { icon: 'github', link: 'https://github.com/emcal/ilreco' },
    ],
  },
})
