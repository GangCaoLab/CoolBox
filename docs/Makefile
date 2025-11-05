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
	@echo ""
	@echo "Building:"
	@echo "  html              build HTML documentation (same as GitHub Actions)"
	@echo "  clean             remove build artifacts"
	@echo "  livehtml          auto-rebuild and serve docs with live reload"
	@echo ""
	@echo "Notebooks:"
	@echo "  list-nb           list all Jupyter notebooks"
	@echo "  update-nb         execute and update all notebooks (sequential)"
	@echo "  update-nb-fast    execute and update all notebooks (parallel)"
	@echo "  update-nb-filter  update notebooks matching FILTER (e.g., FILTER=quick_start)"
	@echo ""
	@echo "Serving:"
	@echo "  serve             serve the built documentation on http://localhost:$(PORT)"
	@echo ""
	@echo "Testing:"
	@echo "  linkcheck         check all external links for integrity"
	@echo "  doctest           run all doctests embedded in the documentation"

.PHONY: help html clean livehtml serve update-nb update-nb-fast update-nb-filter list-nb linkcheck doctest Makefile

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

# Update Jupyter notebooks by executing them
update-nb:
	@echo "Updating Jupyter notebooks (sequential execution)..."
	@python update_notebooks.py --source-dir $(SOURCEDIR)
	@echo ""
	@echo "Notebooks updated successfully!"

# Update notebooks in parallel (faster but uses more resources)
update-nb-fast:
	@echo "Updating Jupyter notebooks (parallel execution)..."
	@python update_notebooks.py --source-dir $(SOURCEDIR) --parallel 4
	@echo ""
	@echo "Notebooks updated successfully!"

# Update specific notebook(s) by name filter
update-nb-filter:
	@if [ -z "$(FILTER)" ]; then \
		echo "Error: Please specify FILTER variable"; \
		echo "Example: make update-nb-filter FILTER=quick_start"; \
		exit 1; \
	fi
	@echo "Updating notebooks matching '$(FILTER)'..."
	@python update_notebooks.py --source-dir $(SOURCEDIR) --filter "$(FILTER)"

# Dry run - list notebooks without executing
list-nb:
	@python update_notebooks.py --source-dir $(SOURCEDIR) --dry-run

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
