# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?= --keep-going
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
PORT          ?= 8000

# Put it first so that "make" without argument is like "make help".
help:
	@echo "Please use 'make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files (same as GitHub Actions)"
	@echo "  clean      to remove build artifacts"
	@echo "  livehtml   to auto-rebuild and serve docs with live reload"
	@echo "  serve      to serve the built documentation"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation"

.PHONY: help html clean livehtml serve linkcheck doctest Makefile

# Build HTML documentation (same as GitHub Actions)
html:
	@echo "Building HTML documentation..."
	rm -rf $(BUILDDIR)/html
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS)
	@echo ""
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
	@echo "Run 'make serve' to view them locally."

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILDDIR)
	@echo "Clean complete."

# Live reload documentation (auto-rebuild on file changes)
livehtml:
	@echo "Starting live documentation server on http://127.0.0.1:$(PORT)"
	@echo "Press Ctrl+C to stop..."
	@command -v sphinx-autobuild >/dev/null 2>&1 || { \
		echo "Error: sphinx-autobuild is not installed."; \
		echo "Install it with: pip install sphinx-autobuild"; \
		exit 1; \
	}
	sphinx-autobuild $(SOURCEDIR) $(BUILDDIR)/html \
		--port $(PORT) \
		--open-browser \
		--watch $(SOURCEDIR) \
		--ignore "*.swp" \
		--ignore "*.swo" \
		--ignore "*~" \
		$(SPHINXOPTS)

# Serve the built documentation
serve:
	@echo "Serving documentation on http://127.0.0.1:$(PORT)"
	@echo "Press Ctrl+C to stop..."
	@cd $(BUILDDIR)/html && python -m http.server $(PORT)

# Check links
linkcheck:
	$(SPHINXBUILD) -b linkcheck $(SOURCEDIR) $(BUILDDIR)/linkcheck $(SPHINXOPTS)
	@echo ""
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

# Run doctests
doctest:
	$(SPHINXBUILD) -b doctest $(SOURCEDIR) $(BUILDDIR)/doctest $(SPHINXOPTS)
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
