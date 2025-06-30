from pathlib import Path

from quantmsio.core.quantms.psm import Psm

TEST_DATA_ROOT = Path(__file__).parent / "examples"


def test_convert_mztab_to_feature():
    mztab_path = TEST_DATA_ROOT / "quantms/dda-lfq-small/PXD040438.mzTab"
    psm = Psm(mztab_path)
    for _ in psm.generate_report():
        print("ok")
