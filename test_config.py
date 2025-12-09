#!/usr/bin/env python3
"""
Test script to validate configuration loading and YAML generation.
This script tests the Boltz-2 SAR Iterator without running actual predictions.
"""

import sys
from pathlib import Path
from boltz2_sar_iterator import (
    load_config_file,
    parse_pocket_residues,
    parse_contact_residues,
    Boltz2SARIterator,
    SARData
)


def test_config_loading():
    """Test loading configuration from JSON and YAML files."""
    print("="*60)
    print("Testing Configuration Loading")
    print("="*60)

    # Test JSON config
    try:
        json_config = load_config_file("config_example.json")
        print("\n✓ Successfully loaded JSON config")
        print(f"  Protein chain: {json_config.get('protein_chain_id')}")
        print(f"  Ligand chain: {json_config.get('ligand_chain_id')}")
        print(f"  Pocket residues: {json_config.get('pocket_residues')}")
        print(f"  Contact residues: {json_config.get('contact_residues')}")
    except Exception as e:
        print(f"\n✗ Failed to load JSON config: {e}")
        return False

    # Test YAML config
    try:
        yaml_config = load_config_file("config_example.yaml")
        print("\n✓ Successfully loaded YAML config")
        print(f"  Protein chain: {yaml_config.get('protein_chain_id')}")
        print(f"  Ligand chain: {yaml_config.get('ligand_chain_id')}")
        print(f"  Pocket residues: {yaml_config.get('pocket_residues')}")
        print(f"  Contact residues: {yaml_config.get('contact_residues')}")
    except Exception as e:
        print(f"\n✗ Failed to load YAML config: {e}")
        return False

    return True


def test_parsing():
    """Test parsing of pocket and contact residues from strings."""
    print("\n" + "="*60)
    print("Testing String Parsing")
    print("="*60)

    # Test pocket residue parsing
    pocket_str = "A:107,A:98,B:50"
    try:
        pocket_residues = parse_pocket_residues(pocket_str)
        print(f"\n✓ Parsed pocket residues: {pocket_str}")
        print(f"  Result: {pocket_residues}")
        assert pocket_residues == [('A', 107), ('A', 98), ('B', 50)]
    except Exception as e:
        print(f"\n✗ Failed to parse pocket residues: {e}")
        return False

    # Test contact residue parsing
    contact_str = "A:20-B:27:5.0,A:21-B:31:5.5"
    try:
        contact_residues = parse_contact_residues(contact_str)
        print(f"\n✓ Parsed contact residues: {contact_str}")
        print(f"  Result: {contact_residues}")
        assert len(contact_residues) == 2
    except Exception as e:
        print(f"\n✗ Failed to parse contact residues: {e}")
        return False

    return True


def test_yaml_generation():
    """Test YAML input file generation."""
    print("\n" + "="*60)
    print("Testing YAML Generation")
    print("="*60)

    # Create a test iterator (without running predictions)
    try:
        iterator = Boltz2SARIterator(
            protein_sequence="MVTPEGNVSLVDESLLVGVTDEDRAV",
            csv_path="example_compounds.csv",
            output_dir="./test_yaml_output",
            protein_chain_id="A",
            ligand_chain_id="L",
            pocket_residues=[("A", 107), ("A", 98)],
            contact_residues=[(("A", 20), ("B", 27), 5.0)],
            template_file=None,
            use_msa_server=False,
            target_r2=0.7,
            max_iterations=1,
            log_level="ERROR"  # Suppress logs for test
        )
        print("\n✓ Created Boltz2SARIterator instance")

        # Create a test compound
        test_compound = SARData(
            compound_id="test_compound",
            smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            activity=6.2
        )

        # Generate YAML
        yaml_path = iterator.create_yaml_input(test_compound, iteration=1)
        print(f"✓ Generated YAML file: {yaml_path}")

        # Read and display the generated YAML
        if yaml_path.exists():
            print("\nGenerated YAML content:")
            print("-" * 60)
            with open(yaml_path, 'r') as f:
                content = f.read()
                print(content)
            print("-" * 60)

            # Validate YAML structure
            import yaml
            with open(yaml_path, 'r') as f:
                yaml_data = yaml.safe_load(f)

            # Check required fields
            assert 'version' in yaml_data
            assert 'sequences' in yaml_data
            assert len(yaml_data['sequences']) == 2
            assert 'protein' in yaml_data['sequences'][0]
            assert 'ligand' in yaml_data['sequences'][1]
            assert yaml_data['sequences'][0]['protein']['id'] == 'A'
            assert yaml_data['sequences'][1]['ligand']['id'] == 'L'
            assert 'constraints' in yaml_data
            assert len(yaml_data['constraints']) == 2  # 1 pocket + 1 contact

            print("✓ YAML structure validation passed")

            # Cleanup
            import shutil
            shutil.rmtree("./test_yaml_output")
            print("✓ Cleaned up test files")

        else:
            print("✗ YAML file was not created")
            return False

    except Exception as e:
        print(f"\n✗ Failed to generate YAML: {e}")
        import traceback
        traceback.print_exc()
        return False

    return True


def test_csv_loading():
    """Test CSV data loading."""
    print("\n" + "="*60)
    print("Testing CSV Loading")
    print("="*60)

    try:
        iterator = Boltz2SARIterator(
            protein_sequence="MVTPEGNVSLVDESLLVGVTDEDRAV",
            csv_path="example_compounds.csv",
            output_dir="./test_csv_output",
            use_msa_server=False,
            target_r2=0.7,
            max_iterations=1,
            log_level="ERROR"
        )

        iterator.load_sar_data()
        print(f"\n✓ Loaded {len(iterator.sar_data)} compounds from CSV")

        # Display first few compounds
        for i, compound in enumerate(iterator.sar_data[:3], 1):
            print(f"\nCompound {i}:")
            print(f"  ID: {compound.compound_id}")
            print(f"  SMILES: {compound.smiles[:50]}...")
            print(f"  Activity: {compound.activity}")

        # Cleanup
        import shutil
        if Path("./test_csv_output").exists():
            shutil.rmtree("./test_csv_output")

    except Exception as e:
        print(f"\n✗ Failed to load CSV: {e}")
        import traceback
        traceback.print_exc()
        return False

    return True


def main():
    """Run all tests."""
    print("\nBoltz-2 SAR Iterator Configuration Test Suite")
    print("="*60)

    tests = [
        ("Config Loading", test_config_loading),
        ("String Parsing", test_parsing),
        ("CSV Loading", test_csv_loading),
        ("YAML Generation", test_yaml_generation),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ Test '{test_name}' crashed: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))

    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)

    all_passed = True
    for test_name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{status}: {test_name}")
        if not passed:
            all_passed = False

    print("\n" + "="*60)
    if all_passed:
        print("All tests PASSED! ✓")
        print("The tool is ready to use.")
        return 0
    else:
        print("Some tests FAILED! ✗")
        print("Please check the errors above.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
