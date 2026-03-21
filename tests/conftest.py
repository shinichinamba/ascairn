# pytest fixtures for ascairn tests
import os
import subprocess
import pytest

# Limit Polars threads to avoid memory issues during tests
os.environ["POLARS_MAX_THREADS"] = "1"

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

# Resource version: set via ASCAIRN_RESOURCE_VERSION env var (default: ver_2024-12-06)
RESOURCE_VERSION = os.environ.get("ASCAIRN_RESOURCE_VERSION", "ver_2024-12-06")


@pytest.fixture(scope="session")
def resource_version():
    return RESOURCE_VERSION


@pytest.fixture(scope="session")
def _ascairn_resource_path():
    """Clone ascairn_resource (devel branch) if not exists and return base path."""
    resource_path = os.path.join(TESTS_DIR, "ascairn_resource")
    if not os.path.exists(resource_path):
        subprocess.run(
            ["git", "clone", "-b", "devel", "https://github.com/friend1ws/ascairn_resource.git"],
            cwd=TESTS_DIR,
            check=True,
        )
    return resource_path


@pytest.fixture(scope="session")
def resource_dir(_ascairn_resource_path):
    """Legacy resource directory."""
    return os.path.join(_ascairn_resource_path, "resource", "legacy", RESOURCE_VERSION)


@pytest.fixture(scope="session")
def panel_dir(_ascairn_resource_path):
    """Panel resource directory (kmer_info/ + hap_info/)."""
    return os.path.join(_ascairn_resource_path, "resource", "panel", "ascairn_paper_2025")


@pytest.fixture(scope="session")
def output_dir():
    path = os.path.join(TESTS_DIR, "output")
    os.makedirs(path, exist_ok=True)
    return path


@pytest.fixture
def expected_dir():
    return os.path.join(TESTS_DIR, "expected", RESOURCE_VERSION)


@pytest.fixture
def s3_cram():
    return "s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram"
