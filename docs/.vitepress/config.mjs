import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'ilreco',
  description: 'Island-clustering reconstruction for square-cell calorimeters',
  base: '/ilreco/',
  themeConfig: {
    nav: [
      { text: 'Install', link: '/install' },
      { text: 'Getting started', link: '/getting-started' },
      { text: 'Algorithms', link: '/algorithms' },
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
        text: 'Algorithms',
        items: [
          { text: 'How ilreco reconstructs showers', link: '/algorithms' },
          { text: 'Clustering & peak finding', link: '/clustering' },
          { text: 'Shower profiles', link: '/profiles' },
          { text: 'The χ² position fit', link: '/chi2-method' },
          { text: 'The log-likelihood method', link: '/log-likelihood' },
          { text: 'Overlapping showers', link: '/two-gamma' },
        ],
      },
      {
        text: 'Development',
        items: [
          { text: 'Building & contributing', link: '/development' },
          { text: 'Generating shower profiles', link: '/generating-profiles' },
          { text: 'Generating likelihood tables', link: '/generating-ll-tables' },
        ],
      },
    ],
    socialLinks: [
      { icon: 'github', link: 'https://github.com/emcal/ilreco' },
    ],
  },
})
