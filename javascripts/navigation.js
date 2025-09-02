// MathJax configuration for PRISM documentation
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

// Custom PRISM documentation scripts
document$.subscribe(function() {
  // ============================================
  // Navigation Auto-Scroll Feature
  // ============================================
  
  // Function to keep active navigation item in view
  function keepActiveNavInView() {
    const sidebar = document.querySelector('.md-sidebar--primary .md-sidebar__scrollwrap');
    const activeItem = document.querySelector('.md-nav--primary .md-nav__link--active');
    
    if (!sidebar || !activeItem) return;
    
    // Get positions
    const sidebarRect = sidebar.getBoundingClientRect();
    const activeRect = activeItem.getBoundingClientRect();
    
    // Calculate if active item is out of view
    const isAbove = activeRect.top < sidebarRect.top;
    const isBelow = activeRect.bottom > sidebarRect.bottom;
    
    // Add indicator classes
    const sidebarElement = document.querySelector('.md-sidebar--primary');
    if (sidebarElement) {
      sidebarElement.classList.toggle('has-active-above', isAbove);
      sidebarElement.classList.toggle('has-active-below', isBelow);
    }
    
    // Auto-scroll if needed
    if (isAbove || isBelow) {
      const offset = 100; // Offset from top/bottom for better visibility
      
      // Calculate scroll position
      const scrollTop = sidebar.scrollTop;
      const itemTop = activeItem.offsetTop;
      const itemHeight = activeItem.offsetHeight;
      const sidebarHeight = sidebar.offsetHeight;
      
      let newScrollTop;
      
      if (isAbove) {
        // Scroll up to show the item with offset from top
        newScrollTop = itemTop - offset;
      } else if (isBelow) {
        // Scroll down to show the item with offset from bottom
        newScrollTop = itemTop - sidebarHeight + itemHeight + offset;
      }
      
      // Smooth scroll to position
      sidebar.scrollTo({
        top: newScrollTop,
        behavior: 'smooth'
      });
      
      // Add animation class
      activeItem.classList.add('auto-scrolled');
      setTimeout(() => {
        activeItem.classList.remove('auto-scrolled');
      }, 600);
    }
  }
  
  // Debounce function for performance
  function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
      const later = () => {
        clearTimeout(timeout);
        func(...args);
      };
      clearTimeout(timeout);
      timeout = setTimeout(later, wait);
    };
  }
  
  // Setup navigation auto-scroll with Intersection Observer
  function setupNavAutoScroll() {
    const sidebar = document.querySelector('.md-sidebar--primary .md-sidebar__scrollwrap');
    if (!sidebar) return;
    
    // Create intersection observer for active items
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.target.classList.contains('md-nav__link--active') && !entry.isIntersecting) {
            // Active item is out of view, scroll it into view
            entry.target.scrollIntoView({
              behavior: 'smooth',
              block: 'center'
            });
          }
        });
      },
      {
        root: sidebar,
        rootMargin: '-50px 0px -50px 0px',
        threshold: 1.0
      }
    );
    
    // Observe all nav links
    const navLinks = document.querySelectorAll('.md-nav--primary .md-nav__link');
    navLinks.forEach(link => observer.observe(link));
  }
  
  // Initialize navigation auto-scroll
  setTimeout(setupNavAutoScroll, 500);
  
  // Handle navigation clicks
  document.querySelectorAll('.md-nav__link').forEach(link => {
    link.addEventListener('click', function() {
      setTimeout(keepActiveNavInView, 300);
    });
  });
  
  // MutationObserver for dynamic content changes
  const navObserver = new MutationObserver(debounce(() => {
    keepActiveNavInView();
  }, 100));
  
  const navElement = document.querySelector('.md-nav--primary');
  if (navElement) {
    navObserver.observe(navElement, {
      attributes: true,
      attributeFilter: ['class'],
      subtree: true
    });
  }
  
  // Handle window resize
  window.addEventListener('resize', debounce(keepActiveNavInView, 200));
  
  // ============================================
  // Code Block Enhancements
  // ============================================
  
  // Add copy functionality to code blocks
  const codeBlocks = document.querySelectorAll('pre > code');
  codeBlocks.forEach(block => {
    // Add line numbers for long code blocks
    const lines = block.textContent.split('\n');
    if (lines.length > 5) {
      block.classList.add('line-numbers');
    }
  });

  // ============================================
  // Smooth Scrolling
  // ============================================
  
  // Smooth scroll for anchor links
  document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
      e.preventDefault();
      const target = document.querySelector(this.getAttribute('href'));
      if (target) {
        target.scrollIntoView({
          behavior: 'smooth',
          block: 'start'
        });
        // Update navigation after scroll
        setTimeout(keepActiveNavInView, 500);
      }
    });
  });

  // ============================================
  // Interactive Features
  // ============================================
  
  // Add interactive tooltips
  const addTooltip = (element, text) => {
    element.setAttribute('data-tooltip', text);
    element.classList.add('has-tooltip');
  };

  // Add tooltips to abbreviations
  document.querySelectorAll('abbr').forEach(abbr => {
    if (abbr.title) {
      addTooltip(abbr, abbr.title);
    }
  });

  // ============================================
  // Table Enhancements
  // ============================================
  
  // Enhance tables with sorting capability
  const tables = document.querySelectorAll('table');
  tables.forEach(table => {
    // Add class for styling
    table.classList.add('sortable');
    
    // Add sort indicators to headers
    const headers = table.querySelectorAll('th');
    headers.forEach((header, index) => {
      header.style.cursor = 'pointer';
      header.setAttribute('data-sort-index', index);
      
      // Add sort icon
      const sortIcon = document.createElement('span');
      sortIcon.className = 'sort-icon';
      sortIcon.innerHTML = ' ⇅';
      header.appendChild(sortIcon);
      
      // Add click handler
      header.addEventListener('click', () => sortTable(table, index));
    });
  });

  // ============================================
  // Progressive Image Loading
  // ============================================
  
  const images = document.querySelectorAll('img[data-src]');
  const imageObserver = new IntersectionObserver((entries, observer) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        const img = entry.target;
        img.src = img.dataset.src;
        img.removeAttribute('data-src');
        imageObserver.unobserve(img);
      }
    });
  });

  images.forEach(img => imageObserver.observe(img));

  // ============================================
  // Feature Card Animations
  // ============================================
  
  const cards = document.querySelectorAll('.grid.cards > *, .feature-card');
  const cardObserver = new IntersectionObserver((entries) => {
    entries.forEach((entry, index) => {
      if (entry.isIntersecting) {
        entry.target.style.setProperty('--card-index', index);
        entry.target.classList.add('animate-in');
      }
    });
  }, { threshold: 0.1 });

  cards.forEach(card => cardObserver.observe(card));

  // ============================================
  // Enhanced Search Experience
  // ============================================
  
  const searchInput = document.querySelector('.md-search__input');
  if (searchInput) {
    // Add search suggestions
    searchInput.addEventListener('input', function() {
      const query = this.value.toLowerCase();
      // Implement search suggestions logic here if needed
    });
  }

  // ============================================
  // Keyboard Shortcuts
  // ============================================
  
  document.addEventListener('keydown', function(e) {
    // Ctrl/Cmd + K for search
    if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
      e.preventDefault();
      const searchInput = document.querySelector('.md-search__input');
      if (searchInput) {
        searchInput.focus();
      }
    }
    
    // Escape to close search
    if (e.key === 'Escape') {
      const searchInput = document.querySelector('.md-search__input');
      if (searchInput && searchInput === document.activeElement) {
        searchInput.blur();
      }
    }
    
    // Arrow keys for navigation
    if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
      const activeLink = document.querySelector('.md-nav__link--active');
      if (activeLink && !e.target.matches('input, textarea')) {
        e.preventDefault();
        const links = Array.from(document.querySelectorAll('.md-nav--primary .md-nav__link'));
        const currentIndex = links.indexOf(activeLink);
        let newIndex;
        
        if (e.key === 'ArrowUp') {
          newIndex = Math.max(0, currentIndex - 1);
        } else {
          newIndex = Math.min(links.length - 1, currentIndex + 1);
        }
        
        if (links[newIndex]) {
          links[newIndex].click();
        }
      }
    }
  });

  // ============================================
  // Reading Progress Indicator
  // ============================================
  
  const article = document.querySelector('article');
  if (article) {
    const progressBar = document.createElement('div');
    progressBar.className = 'reading-progress';
    progressBar.style.cssText = `
      position: fixed;
      top: 0;
      left: 0;
      height: 3px;
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      z-index: 1000;
      transition: width 0.2s ease;
    `;
    document.body.appendChild(progressBar);

    window.addEventListener('scroll', () => {
      const scrollTop = window.scrollY;
      const scrollHeight = article.scrollHeight - window.innerHeight;
      const progress = (scrollTop / scrollHeight) * 100;
      progressBar.style.width = `${progress}%`;
    });
  }

  // ============================================
  // Code Syntax Highlighting
  // ============================================
  
  if (typeof Prism !== 'undefined') {
    Prism.highlightAll();
  }

  // ============================================
  // Mermaid Diagrams
  // ============================================
  
  if (typeof mermaid !== 'undefined') {
    mermaid.initialize({
      startOnLoad: true,
      theme: 'default',
      themeVariables: {
        primaryColor: '#4051b5',
        primaryTextColor: '#fff',
        primaryBorderColor: '#7c4dff',
        lineColor: '#667eea',
        secondaryColor: '#764ba2',
        tertiaryColor: '#9575cd'
      }
    });
  }
});

// ============================================
// Table Sorting Function
// ============================================

function sortTable(table, columnIndex) {
  const tbody = table.querySelector('tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const header = table.querySelectorAll('th')[columnIndex];
  const isAscending = header.classList.contains('sort-asc');
  
  // Sort rows
  rows.sort((a, b) => {
    const aValue = a.cells[columnIndex].textContent.trim();
    const bValue = b.cells[columnIndex].textContent.trim();
    
    // Try to parse as numbers
    const aNum = parseFloat(aValue);
    const bNum = parseFloat(bValue);
    
    if (!isNaN(aNum) && !isNaN(bNum)) {
      return isAscending ? bNum - aNum : aNum - bNum;
    }
    
    // Sort as strings
    return isAscending 
      ? bValue.localeCompare(aValue)
      : aValue.localeCompare(bValue);
  });
  
  // Update table
  tbody.innerHTML = '';
  rows.forEach(row => tbody.appendChild(row));
  
  // Update sort indicators
  table.querySelectorAll('th').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
  });
  
  header.classList.add(isAscending ? 'sort-desc' : 'sort-asc');
}

// ============================================
// Custom Animations Style Injection
// ============================================

const style = document.createElement('style');
style.textContent = `
  @keyframes fadeInUp {
    from {
      opacity: 0;
      transform: translateY(20px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
  
  .animate-in {
    animation: fadeInUp 0.5s ease forwards;
  }
  
  .has-tooltip {
    position: relative;
    cursor: help;
  }
  
  .has-tooltip:hover::after {
    content: attr(data-tooltip);
    position: absolute;
    bottom: 100%;
    left: 50%;
    transform: translateX(-50%);
    padding: 0.5rem;
    background: rgba(0, 0, 0, 0.8);
    color: white;
    border-radius: 4px;
    font-size: 0.875rem;
    white-space: nowrap;
    z-index: 1000;
    pointer-events: none;
  }
  
  .sort-asc .sort-icon::after {
    content: ' ↑';
  }
  
  .sort-desc .sort-icon::after {
    content: ' ↓';
  }
`;
document.head.appendChild(style);