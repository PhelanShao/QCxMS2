#!/usr/bin/env python3
"""
Test script to verify all QCxMS2 Python modules can be imported correctly.
"""

import sys
from pathlib import Path

# Add src directory to Python path
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

def test_imports():
    """Test importing all QCxMS2 modules"""
    print("Testing QCxMS2 Python module imports...")
    
    try:
        print("  Importing constants...", end=" ")
        from qcxms import constants
        print("✓")
        
        print("  Importing data...", end=" ")
        from qcxms import data
        print("✓")
        
        print("  Importing iomod...", end=" ")
        from qcxms import iomod
        print("✓")
        
        print("  Importing utility...", end=" ")
        from qcxms import utility
        print("✓")
        
        print("  Importing argparser...", end=" ")
        from qcxms import argparser
        print("✓")
        
        print("  Importing qmmod...", end=" ")
        from qcxms import qmmod
        print("✓")
        
        print("  Importing charges...", end=" ")
        from qcxms import charges
        print("✓")
        
        print("  Importing fragmentation...", end=" ")
        from qcxms import fragmentation
        print("✓")
        
        print("  Importing mcsimu...", end=" ")
        from qcxms import mcsimu
        print("✓")
        
        print("  Importing reaction...", end=" ")
        from qcxms import reaction
        print("✓")
        
        print("  Importing tsmod...", end=" ")
        from qcxms import tsmod
        print("✓")
        
        print("  Importing cid...", end=" ")
        from qcxms import cid
        print("✓")
        
        print("  Importing plotting...", end=" ")
        from qcxms import plotting
        print("✓")
        
        print("  Importing main...", end=" ")
        from qcxms import main
        print("✓")
        
        print("\n✅ All modules imported successfully!")
        return True
        
    except ImportError as e:
        print(f"\n❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        return False

def test_basic_functionality():
    """Test basic functionality of key modules"""
    print("\nTesting basic functionality...")
    
    try:
        from qcxms import constants, utility, data
        
        # Test constants
        print(f"  Constants - AUTOEV: {constants.AUTOEV}")
        print(f"  Constants - KB_EV_K: {constants.KB_EV_K}")
        
        # Test utility functions
        temp = utility.calctemp(30, 1.5)
        print(f"  Utility - calctemp(30, 1.5): {temp:.2f} K")
        
        # Test data structures
        env = data.RunTypeData()
        print(f"  Data - Default charge: {env.chrg}")
        print(f"  Data - Default mode: {env.mode}")
        
        print("✅ Basic functionality tests passed!")
        return True
        
    except Exception as e:
        print(f"❌ Functionality test error: {e}")
        return False

def main():
    """Main test function"""
    print("QCxMS2 Python Implementation - Module Test")
    print("=" * 50)
    
    success = True
    
    # Test imports
    if not test_imports():
        success = False
    
    # Test basic functionality
    if not test_basic_functionality():
        success = False
    
    print("\n" + "=" * 50)
    if success:
        print("🎉 All tests passed! QCxMS2 Python implementation is ready.")
        return 0
    else:
        print("💥 Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())