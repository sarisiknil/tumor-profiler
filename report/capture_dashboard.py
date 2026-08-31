#!/usr/bin/env python3
"""Start the dashboard, screenshot each page at presentation resolution, and shut it down.

Screenshots are written to report/figures/. They are taken from the live application, so what appears in the
report and the presentation is what the application actually renders — not a mock-up.

  python3 report/capture_dashboard.py --results results
"""
import argparse, socket, subprocess, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PAGES = [("", "01_home"), ("QC_and_Library", "02_qc_library"), ("Variants", "03_variants"),
         ("RNA_and_Fusions", "04_rna_fusions"), ("Pathways", "05_pathways"),
         ("Therapies", "06_therapies"), ("Method_and_Limits", "07_method")]


def free_port():
    s = socket.socket()
    s.bind(("", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--outdir", default="report/figures")
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=1000)
    ap.add_argument("--scale", type=int, default=2, help="device pixel ratio; 2 gives retina-quality images")
    a = ap.parse_args()
    out = REPO / a.outdir
    out.mkdir(parents=True, exist_ok=True)
    port = free_port()

    proc = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "dashboard/app.py",
         "--server.port", str(port), "--server.headless", "true",
         "--browser.gatherUsageStats", "false", "--theme.base", "dark",
         "--", "--results", a.results],
        cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(60):
            try:
                socket.create_connection(("127.0.0.1", port), timeout=1).close()
                break
            except OSError:
                time.sleep(1)
        time.sleep(6)

        from playwright.sync_api import sync_playwright
        with sync_playwright() as pw:
            # Use the Chrome already installed on the machine rather than downloading a private copy:
            # the disk is nearly full, and a second browser build is 500 MB for no benefit.
            try:
                browser = pw.chromium.launch(channel="chrome")
            except Exception:
                browser = pw.chromium.launch()
            page = browser.new_page(viewport={"width": a.width, "height": a.height},
                                    device_scale_factor=a.scale)
            for path, name in PAGES:
                page.goto(f"http://127.0.0.1:{port}/{path}", wait_until="networkidle", timeout=90000)
                page.wait_for_timeout(4500)
                # hide Streamlit's own chrome so the figure shows the content, not the framework
                page.add_style_tag(content="""
                    header[data-testid='stHeader'], #MainMenu, footer,
                    [data-testid='stToolbar'], [data-testid='stDecoration'],
                    [data-testid='stSidebarCollapsedControl'] { display: none !important; }
                    section[data-testid='stMain'] { padding-top: 0 !important; }
                """)
                page.wait_for_timeout(1200)
                page.screenshot(path=str(out / f"{name}.png"), full_page=False)
                page.screenshot(path=str(out / f"{name}_full.png"), full_page=True)
                print(f"  captured {name}")
            browser.close()
    finally:
        proc.terminate()
        proc.wait(timeout=20)
    print(f"\nwrote {len(list(out.glob('*.png')))} images to {out}")


if __name__ == "__main__":
    main()
