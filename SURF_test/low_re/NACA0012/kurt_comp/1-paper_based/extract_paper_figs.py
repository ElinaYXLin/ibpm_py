"""
extract_paper_figs.py

Crops the paper's own wake-vorticity figures (Kurtulus 2019, Figure 2 =
steady, Figure 6 = f=4Hz) at alpha_0 = 0/9/12 deg -- the same three angles
gen_kurt_figs.py's wake_contours() plots for py_static/cpp_static -- so
wake_steady.png / wake_f4hz.png can show what's actually being compared
against, not just ibpm's own output. Renders the PDF pages with PyMuPDF
(poppler/pdftoppm isn't installed on this machine) and crops each row by
pixel offsets found once by inspection of the rendered page (hardcoded
below -- brittle if the PDF is ever swapped for a different version, but
this is a one-time extraction, not something run as part of the regular
pipeline).

Usage: python3 extract_paper_figs.py
Output: paper_figs/{steady,f4hz}_a{00,09,12}.png
"""
import pathlib
import fitz

PDF = pathlib.Path(
    "/Users/elina/Downloads/"
    "kurtulus-2019-unsteady-aerodynamics-of-a-pitching-naca-0012-airfoil-at-low-reynolds-number.pdf"
)
OUT = pathlib.Path(__file__).resolve().parent / "paper_figs"
OUT.mkdir(exist_ok=True)

# (pdf page index (0-based), crop box in 200dpi-rendered pixel coords)
# Figure 2 (page 4, steady, instantaneous vorticity): rows for
# alpha_0 = 0,7,8,9,10,11,12 deg -- only 0/9/12 kept, matching
# gen_kurt_figs.py's wake_contours() angle set.
STEADY_CROPS = {
    "a00": (3, (340, 860, 1300, 990)),
    "a09": (3, (340, 1157, 1300, 1243)),
    "a12": (3, (340, 1409, 1300, 1497)),
}
# Figure 6 (page 6, f=4Hz, instantaneous vorticity, t=36s): rows for
# alpha_0 = 0,1,2,3,9,10,11,12 deg -- only 0/9/12 kept.
F4HZ_CROPS = {
    "a00": (5, (340, 200, 1300, 280)),
    "a09": (5, (340, 512, 1300, 610)),
    "a12": (5, (340, 786, 1300, 882)),
}


def main():
    doc = fitz.open(PDF)
    page_cache = {}
    for label, crops in (("steady", STEADY_CROPS), ("f4hz", F4HZ_CROPS)):
        for tag, (page_idx, box) in crops.items():
            if page_idx not in page_cache:
                page_cache[page_idx] = doc[page_idx].get_pixmap(dpi=200)
            pix = page_cache[page_idx]
            # crop via PIL since fitz.Pixmap has no direct crop
            import io
            from PIL import Image
            im = Image.open(io.BytesIO(pix.tobytes("png")))
            out_path = OUT / f"{label}_{tag}.png"
            im.crop(box).save(out_path)
            print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
