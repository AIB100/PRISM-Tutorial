window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    tags: 'ams',
    packages: {
      '[+]': ['mathtools', 'cases', 'braket', 'mhchem']
    }
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  },
  loader: {
    load: ['[tex]/mathtools', '[tex]/cases', '[tex]/braket', '[tex]/mhchem']
  },
  chtml: {
    scale: 1.0,
    displayAlign: 'left',
    displayIndent: '2em'
  },
  svg: {
    scale: 1.0,
    displayAlign: 'left',
    displayIndent: '2em'
  }
};
