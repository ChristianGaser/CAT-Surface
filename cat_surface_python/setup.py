"""
Setup script for CAT-Surface Python wrapper.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cat-surface-python",
    version="0.1.0",
    author="GitHub Copilot",
    author_email="",
    description="Python wrapper for CAT-Surface neuroimaging tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ChristianGaser/CAT-Surface",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        # No external dependencies - uses only Python standard library
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
    },
    entry_points={
        "console_scripts": [
            # Could add command-line interfaces here
        ],
    },
    keywords="neuroimaging surface analysis cortex brain CAT-Surface",
    project_urls={
        "Bug Reports": "https://github.com/ChristianGaser/CAT-Surface/issues",
        "Source": "https://github.com/ChristianGaser/CAT-Surface",
        "Documentation": "https://github.com/ChristianGaser/CAT-Surface",
    },
)