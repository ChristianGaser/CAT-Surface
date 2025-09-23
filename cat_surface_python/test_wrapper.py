"""
Simple test script for CAT-Surface Python wrapper.

This script tests basic functionality without requiring actual data files.
"""

import sys
import os

# Add the cat_surface module to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import cat_surface
from cat_surface import CATConfig, CATSurfaceError, CATExecutableNotFoundError


def test_configuration():
    """Test configuration and executable detection."""
    print("Testing configuration...")
    
    try:
        # Test auto-detection
        config = CATConfig()
        print(f"✓ Auto-detected executable path: {config.executable_path}")
        
        # Test a few executables
        test_tools = ["CAT_SurfConvert", "CAT_VolAmap", "CAT_SurfDistance"]
        for tool in test_tools:
            try:
                exe_path = config.get_executable(tool)
                print(f"✓ Found {tool}: {exe_path}")
            except CATExecutableNotFoundError:
                print(f"✗ Missing {tool}")
        
        return True
        
    except CATExecutableNotFoundError as e:
        print(f"✗ Configuration failed: {e}")
        print("Please ensure CAT-Surface is compiled and available.")
        return False


def test_error_handling():
    """Test error handling with invalid inputs."""
    print("\nTesting error handling...")
    
    try:
        # This should fail gracefully
        result = cat_surface.segment_with_amap(
            input_file="/nonexistent/file.nii",
            label_file="/nonexistent/labels.nii"
        )
        print("✗ Expected error but got success")
        return False
        
    except FileNotFoundError:
        print("✓ Correctly caught FileNotFoundError for missing input")
        return True
        
    except Exception as e:
        print(f"✗ Unexpected error type: {type(e).__name__}: {e}")
        return False


def test_imports():
    """Test that all modules import correctly."""
    print("\nTesting imports...")
    
    try:
        from cat_surface.volume import VolAmap, VolThicknessPbt
        from cat_surface.surface import SurfAverage, SurfDistance
        print("✓ All modules imported successfully")
        
        # Test class instantiation
        vol_tool = VolAmap()
        surf_tool = SurfAverage()
        print("✓ Tool classes instantiated successfully")
        
        return True
        
    except Exception as e:
        print(f"✗ Import error: {e}")
        return False


def test_help_functionality():
    """Test that tools can show help/usage without crashing."""
    print("\nTesting help functionality...")
    
    try:
        config = cat_surface.get_default_config()
        
        # Try to run a tool with no arguments (should show usage)
        from cat_surface.volume import VolAmap
        tool = VolAmap(config)
        
        # This will fail but should give us info about the tool
        try:
            result = tool._run_command([])  # No arguments
        except CATSurfaceError as e:
            print(f"✓ Tool responded appropriately to no arguments: {type(e).__name__}")
            return True
        
        print("✓ Tool executed without error (unexpected but OK)")
        return True
        
    except CATExecutableNotFoundError:
        print("? Cannot test help - executables not found")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def main():
    """Run all tests."""
    print("CAT-Surface Python Wrapper Test Suite")
    print("=" * 40)
    
    tests = [
        test_imports,
        test_configuration,
        test_error_handling,
        test_help_functionality,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"✗ Test {test.__name__} crashed: {e}")
            results.append(False)
    
    print("\n" + "=" * 40)
    print("Test Results:")
    passed = sum(results)
    total = len(results)
    
    for i, (test, result) in enumerate(zip(tests, results)):
        status = "PASS" if result else "FAIL"
        print(f"  {test.__name__}: {status}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The wrapper appears to be working correctly.")
    elif passed > 0:
        print("⚠️  Some tests passed. Check configuration and CAT-Surface installation.")
    else:
        print("❌ All tests failed. Please check installation and configuration.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)