import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

// This runs in Node.js - Don't use client-side code here (browser APIs, JSX...)

const config: Config = {
  title: 'Espresso_LU',
  tagline: 'Notes and tutorials on first-principles calculations with Quantum ESPRESSO',
  favicon: 'img/uni_icon.ico',

  // === 部署相关 ===
  url: 'https://zixiantan.github.io',      // 你的 GitHub Pages 主域
  baseUrl: '/Espresso_LU/',                // 仓库名
  organizationName: 'ZixianTAN',           // GitHub 用户名
  projectName: 'Espresso_LU',              // 仓库名
  deploymentBranch: 'gh-pages',
  trailingSlash: false,

  //onBrokenLinks: 'throw',
  //onBrokenMarkdownLinks: 'warn',

  // Future flags, see https://docusaurus.io/docs/api/docusaurus-config#future
  //future: {
  //  v4: true, // Improve compatibility with the upcoming Docusaurus v4
  //},

  // Set the production url of your site here
  //url: 'https://your-docusaurus-site.example.com',
  // Set the /<baseUrl>/ pathname under which your site is served
  // For GitHub pages deployment, it is often '/<projectName>/'
  //baseUrl: '/',

  // GitHub pages deployment config.
  // If you aren't using GitHub pages, you don't need these.
  //organizationName: 'facebook', // Usually your GitHub org/user name.
  //projectName: 'docusaurus', // Usually your repo name.

  //onBrokenLinks: 'throw',

  // Even if you don't use internationalization, you can use this field to set
  // useful metadata like html lang. For example, if your site is Chinese, you
  // may want to replace "en" with "zh-Hans".
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      {
        docs: {
          path: '../docs',
          routeBasePath: '/', 
          sidebarPath: './sidebars.ts',
          editUrl: 'https://github.com/ZixianTAN/Espresso_LU/edit/main/docs/',
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
          //sidebarPath: './sidebars.ts',
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          //editUrl:
          //  'https://github.com/facebook/docusaurus/tree/main/packages/create-docusaurus/templates/shared/',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    // Replace with your project's social card
    image: 'img/docusaurus-social-card.jpg',
    colorMode: {
      respectPrefersColorScheme: true,
    },
    navbar: {
      title: 'Espresso_LU',
      logo: {
        alt: 'Lancaster University logo',
        src: '/img/uni_logo_full.svg',
      },
      items: [
        {
          type: 'docSidebar',
          sidebarId: 'docs',
          position: 'left',
          label: 'Docs',
        },
        //{to: '/blog', label: 'Blog', position: 'left'},
        {
          href: 'https://github.com/ZixianTAN/Espresso_LU',
          label: 'GitHub',
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
            { label: 'Foundations', to: 'foundations' },
            { label: 'Quantum ESPRESSO', to: 'qe' },
            { label: 'Softwares', to: 'softwares' },
            { label: 'HPC', to: 'hpc' },
          ],
        },
        //{
        //  title: 'Community',
        //  items: [
        //    {
        //      label: 'Stack Overflow',
        //      href: 'https://stackoverflow.com/questions/tagged/docusaurus',
        //    },
        //    {
        //      label: 'Discord',
        //      href: 'https://discordapp.com/invite/docusaurus',
        //    },
        //    {
        //      label: 'X',
        //      href: 'https://x.com/docusaurus',
        //    },
        //  ],
        //},
        {
          title: 'More',
          items: [
            //{
            //  label: 'Blog',
            //  to: '/blog',
            //},
            { label: 'GitHub', href: 'https://github.com/ZixianTAN/Espresso_LU' },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} Zixian Tan.`,
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
