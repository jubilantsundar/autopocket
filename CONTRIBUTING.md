# Contributing to AutoPocket

Thank you for your interest in contributing to AutoPocket! This document provides guidelines for contributing.

## Ways to Contribute

- 🐛 **Report bugs** via [GitHub Issues](https://github.com/YOUR_USERNAME/autopocket/issues)
- 💡 **Suggest features** via [GitHub Issues](https://github.com/YOUR_USERNAME/autopocket/issues)
- 📝 **Improve documentation**
- 🧪 **Add tests**
- 💻 **Submit code**

## Getting Started

### 1. Fork and Clone

```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/autopocket.git
cd autopocket

# Add upstream remote
git remote add upstream https://github.com/ORIGINAL_OWNER/autopocket.git
```

### 2. Set Up Development Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install development dependencies
pip install -e ".[dev]"

# Install Boltz-2
pip install boltz[cuda] -U
```

### 3. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

## Development Guidelines

### Code Style

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Use [Black](https://github.com/psf/black) for formatting: `black .`
- Use [flake8](https://flake8.pycqa.org/) for linting: `flake8 .`
- Add type hints where appropriate
- Write docstrings for all functions/classes

Example:
```python
def calculate_correlation(
    experimental: List[float],
    predicted: List[float]
) -> Dict[str, float]:
    """
    Calculate correlation metrics between experimental and predicted values.

    Args:
        experimental: List of experimental activity values
        predicted: List of predicted affinity values

    Returns:
        Dictionary with correlation metrics (r2, spearman, etc.)

    Raises:
        ValueError: If arrays have different lengths
    """
    pass
```

### Testing

- Add tests for new features
- Ensure all tests pass: `pytest`
- Aim for >80% code coverage

```bash
# Run tests
pytest

# Run with coverage
pytest --cov=autopocket --cov-report=html
```

### Documentation

- Update README.md if adding new features
- Add docstrings to all functions
- Update QUICKSTART.md for user-facing changes
- Add examples to `examples/` directory

## Commit Guidelines

### Commit Messages

Use clear, descriptive commit messages:

```
feat: add ROC-AUC metric calculation
fix: correct directory structure for Boltz-2 v2.2+
docs: update installation instructions for M1 Macs
test: add unit tests for correlation metrics
refactor: simplify YAML generation logic
```

**Format:**
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `test:` - Adding/updating tests
- `refactor:` - Code refactoring
- `perf:` - Performance improvements
- `chore:` - Maintenance tasks

### Making Commits

```bash
git add .
git commit -m "feat: add your feature description"
git push origin feature/your-feature-name
```

## Pull Request Process

### 1. Before Submitting

- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] Documentation is updated
- [ ] Commit messages are clear
- [ ] Branch is up to date with main

```bash
# Update your branch
git fetch upstream
git rebase upstream/main
```

### 2. Submit PR

1. Go to GitHub and create a Pull Request
2. Fill out the PR template
3. Link any related issues
4. Request review

### 3. PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
Describe the tests you ran

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Comments added for complex code
- [ ] Documentation updated
- [ ] No new warnings
- [ ] Tests added/updated
- [ ] All tests passing
```

### 4. Review Process

- Maintainers will review your PR
- Address any requested changes
- Once approved, your PR will be merged!

## Code Review Guidelines

When reviewing others' PRs:

- Be respectful and constructive
- Focus on the code, not the person
- Suggest improvements clearly
- Approve when satisfied

## Community

- Be respectful and inclusive
- Follow the [Code of Conduct](CODE_OF_CONDUCT.md)
- Help others in discussions and issues

## Questions?

- Open a [Discussion](https://github.com/YOUR_USERNAME/autopocket/discussions)
- Email: your.email@example.com

---

**Thank you for contributing to AutoPocket! 🙏**
