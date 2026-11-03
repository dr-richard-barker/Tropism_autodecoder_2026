/**
 * Viewer helpers: PDF fullscreen toggle, Figure Gallery Lightbox, and BibTeX copy
 */
(function() {
  'use strict';

  // 1. PDF Fullscreen Toggle
  const pdfWrap = document.getElementById('pdf-wrap');
  const fullscreenBtn = document.getElementById('pdf-fullscreen-btn');
  const fullscreenExit = document.getElementById('pdf-fullscreen-exit');

  if (fullscreenBtn && pdfWrap) {
    fullscreenBtn.addEventListener('click', () => {
      pdfWrap.classList.add('fullscreen');
      document.body.style.overflow = 'hidden';
    });
  }

  if (fullscreenExit && pdfWrap) {
    fullscreenExit.addEventListener('click', () => {
      pdfWrap.classList.remove('fullscreen');
      document.body.style.overflow = '';
    });
  }

  // 2. Lightbox Modal for Publication Figures
  const modal = document.getElementById('lightbox-modal');
  const modalImg = document.getElementById('lightbox-img');
  const modalTitle = document.getElementById('lightbox-title');
  const modalDesc = document.getElementById('lightbox-desc');
  const modalOpenBtn = document.getElementById('lightbox-open-btn');
  const modalClose = document.getElementById('lightbox-close');

  if (modal && modalImg) {
    document.querySelectorAll('.gallery-card, .diagram-box img').forEach(card => {
      card.addEventListener('click', (e) => {
        const img = card.querySelector('img') || card;
        const title = card.getAttribute('data-title') || img.alt || 'Publication Figure';
        const desc = card.getAttribute('data-desc') || '';
        
        modalImg.src = img.src;
        modalImg.alt = title;
        if (modalTitle) modalTitle.textContent = title;
        if (modalDesc) modalDesc.textContent = desc;
        if (modalOpenBtn) modalOpenBtn.href = img.src;
        modal.classList.add('active');
        document.body.style.overflow = 'hidden';
      });
    });

    const closeModal = () => {
      modal.classList.remove('active');
      document.body.style.overflow = '';
    };

    if (modalClose) modalClose.addEventListener('click', closeModal);
    modal.addEventListener('click', (e) => {
      if (e.target === modal || e.target.classList.contains('lightbox-img-wrap')) {
        closeModal();
      }
    });

    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        if (modal.classList.contains('active')) closeModal();
        if (pdfWrap && pdfWrap.classList.contains('fullscreen')) {
          pdfWrap.classList.remove('fullscreen');
          document.body.style.overflow = '';
        }
      }
    });
  }

  // 3. Copy BibTeX
  const copyBtn = document.getElementById('copy-bibtex-btn');
  const bibtexCode = document.getElementById('bibtex-code');

  if (copyBtn && bibtexCode) {
    copyBtn.addEventListener('click', async () => {
      try {
        await navigator.clipboard.writeText(bibtexCode.textContent.trim());
        const originalText = copyBtn.textContent;
        copyBtn.textContent = '✓ Copied!';
        copyBtn.style.borderColor = 'var(--blue)';
        copyBtn.style.color = 'var(--blue)';
        setTimeout(() => {
          copyBtn.textContent = originalText;
          copyBtn.style.borderColor = '';
          copyBtn.style.color = '';
        }, 2200);
      } catch (err) {
        console.error('Clipboard copy failed:', err);
      }
    });
  }
})();
