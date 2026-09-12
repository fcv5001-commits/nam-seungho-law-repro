"""Generate the public quantitative whitepaper and residual plots."""

from __future__ import annotations

import json
import urllib.request
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, KeepTogether, PageBreak, PageTemplate,
    Paragraph, Spacer, Table, TableStyle,
)

from data_loader import SOURCES, fetch
from fit_pantheon_ns import chi2_profiled, distance_modulus, load_sample

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
ASSETS = ROOT / "whitepaper_assets"
OUTPUT = ROOT / "Nam_Seung_ho_Law_Open_Data_Whitepaper_v1.0.pdf"
FONT_DIR = ROOT / "data" / "fonts"
FONT_PATH = FONT_DIR / "NotoSansKR.ttf"
FONT_URL = "https://raw.githubusercontent.com/google/fonts/main/ofl/notosanskr/NotoSansKR%5Bwght%5D.ttf"


def ensure_font():
    FONT_DIR.mkdir(parents=True, exist_ok=True)
    if not FONT_PATH.exists():
        request = urllib.request.Request(FONT_URL, headers={"User-Agent": "nam-seungho-law-repro/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            FONT_PATH.write_bytes(response.read())
    return FONT_PATH


def load_json(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def make_plots():
    font_path = ensure_font()
    from matplotlib import font_manager
    font_prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams["font.family"] = font_prop.get_name()
    plt.rcParams["axes.unicode_minus"] = False
    ASSETS.mkdir(exist_ok=True)
    pantheon = load_json("pantheon_fit.json")
    data_path, _ = fetch("pantheon_plus", ROOT / "data")
    z, observed, sigma = load_sample(data_path)

    om = pantheon["baseline_no_lambda"]["omega_m"]
    alpha = pantheon["ns_candidate"]["alpha"]
    nu = pantheon["ns_candidate"]["nu"]
    e_base = lambda zz: np.sqrt(1 + om * ((1 + zz) ** 3 - 1))
    e_ns = lambda zz: np.exp(0.5 * alpha * ((1 + zz) ** nu - 1))
    mu_base = distance_modulus(z, e_base) + pantheon["baseline_no_lambda"]["offset"]
    mu_ns = distance_modulus(z, e_ns) + pantheon["ns_candidate"]["offset"]
    res_base = observed - mu_base
    res_ns = observed - mu_ns

    order = np.argsort(z)
    edges = np.quantile(z, np.linspace(0, 1, 21))
    centers, bbase, bns, berr = [], [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (z >= left) & (z <= right)
        centers.append(np.median(z[mask]))
        bbase.append(np.average(res_base[mask], weights=1 / sigma[mask] ** 2))
        bns.append(np.average(res_ns[mask], weights=1 / sigma[mask] ** 2))
        berr.append(np.sqrt(1 / np.sum(1 / sigma[mask] ** 2)))

    plt.figure(figsize=(9, 4.8))
    plt.scatter(z[order], res_base[order], s=5, alpha=.12, color="#78909c", label="individual residuals")
    plt.errorbar(centers, bbase, yerr=berr, fmt="o-", lw=1.4, ms=4, color="#e0a629", label="no-Lambda baseline")
    plt.errorbar(centers, bns, yerr=berr, fmt="o-", lw=1.4, ms=4, color="#00a99d", label="NS candidate")
    plt.axhline(0, color="black", lw=.8)
    plt.xscale("log"); plt.xlabel("redshift zHD"); plt.ylabel("magnitude residual")
    plt.title("Pantheon+ diagonal-error residuals (20 quantile bins)")
    plt.legend(frameon=False, ncol=3, fontsize=8); plt.grid(alpha=.15); plt.tight_layout()
    pantheon_plot = ASSETS / "pantheon_residuals.png"
    plt.savefig(pantheon_plot, dpi=180); plt.close()

    cmb_path, _ = fetch("planck_tt_binned", ROOT / "data")
    table = np.loadtxt(cmb_path, comments="#")
    ell, obs, elo, ehi, best = table.T
    norm = (obs - best) / (0.5 * (elo + ehi))
    plt.figure(figsize=(9, 4.8))
    plt.axhspan(-2, 2, color="#00a99d", alpha=.08, label="±2σ band")
    plt.axhline(0, color="black", lw=.8)
    plt.vlines(ell, 0, norm, color="#486878", lw=.8, alpha=.7)
    plt.scatter(ell, norm, s=16, color="#e0a629")
    plt.xlabel("multipole ell"); plt.ylabel("(observed - BestFit) / sigma")
    plt.title("Planck 2018 TT normalized residuals (83 released bins)")
    plt.grid(alpha=.15); plt.legend(frameon=False); plt.tight_layout()
    cmb_plot = ASSETS / "planck_tt_residuals.png"
    plt.savefig(cmb_plot, dpi=180); plt.close()
    return pantheon_plot, cmb_plot


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#23404d")); canvas.line(18*mm, 14*mm, 192*mm, 14*mm)
    canvas.setFont("NotoKR", 8); canvas.setFillColor(colors.HexColor("#40545e"))
    canvas.drawString(18*mm, 9*mm, "남승호법칙 공개 오픈데이터 백서 v1.0")
    canvas.drawRightString(192*mm, 9*mm, f"{doc.page}")
    canvas.restoreState()


def build_pdf():
    pantheon_plot, cmb_plot = make_plots()
    p = load_json("pantheon_fit.json"); c = load_json("cmb_residual.json")
    pdfmetrics.registerFont(TTFont("NotoKR", str(ensure_font())))
    styles = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=styles["BodyText"], fontName="NotoKR", fontSize=10.5, leading=17, textColor=colors.HexColor("#18313d"), spaceAfter=8)
    small = ParagraphStyle("small", parent=body, fontSize=8.5, leading=13, textColor=colors.HexColor("#49616b"))
    h1 = ParagraphStyle("h1", parent=body, fontSize=25, leading=34, textColor=colors.HexColor("#07121c"), spaceAfter=15)
    h2 = ParagraphStyle("h2", parent=body, fontSize=17, leading=23, textColor=colors.HexColor("#007f78"), spaceBefore=8, spaceAfter=12)
    center = ParagraphStyle("center", parent=body, alignment=TA_CENTER)
    label = ParagraphStyle("label", parent=body, fontSize=9, textColor=colors.HexColor("#007f78"), alignment=TA_CENTER)

    doc = BaseDocTemplate(str(OUTPUT), pagesize=A4, leftMargin=18*mm, rightMargin=18*mm, topMargin=18*mm, bottomMargin=19*mm, title="남승호법칙 공개 오픈데이터 정량검증 백서 v1.0", author="남승호")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="main")
    doc.addPageTemplates(PageTemplate(id="all", frames=frame, onPage=footer))
    story = []
    story += [Spacer(1, 25*mm), Paragraph("NAM SEUNG-HO LAW · OPEN DATA RECORD", label), Spacer(1, 8*mm), Paragraph("남승호법칙<br/>공개 오픈데이터 정량검증 백서", h1), Paragraph("Residual Plots · Reproduction Code · SHA-256", center), Spacer(1, 18*mm)]
    status = Table([["공개자료 계산 재현", "PASS"], ["독립 물리검증", "OPEN · 미봉인"]], colWidths=[88*mm, 58*mm])
    status.setStyle(TableStyle([("FONTNAME",(0,0),(-1,-1),"NotoKR"),("FONTSIZE",(0,0),(-1,-1),11),("BACKGROUND",(0,0),(0,-1),colors.HexColor("#e8f2f1")),("BACKGROUND",(1,0),(1,0),colors.HexColor("#d9f3ed")),("BACKGROUND",(1,1),(1,1),colors.HexColor("#fff1d7")),("GRID",(0,0),(-1,-1),.6,colors.HexColor("#73909a")),("ALIGN",(1,0),(1,-1),"CENTER"),("VALIGN",(0,0),(-1,-1),"MIDDLE"),("TOPPADDING",(0,0),(-1,-1),10),("BOTTOMPADDING",(0,0),(-1,-1),10)]))
    story += [status, Spacer(1, 18*mm), Paragraph("저자·발명자: 남승호", center), Paragraph(f"공개판 v1.0 · {date.today().isoformat()} · 비동료심사 연구기록", small), PageBreak()]

    story += [Paragraph("1. 요약", h2), Paragraph("이 백서는 남승호법칙 우주론 후보식에 대해 공개된 Pantheon+ 초신성 자료와 Planck 2018 TT 구간자료를 사용한 재현 가능한 계산 기록이다. 코드, 입력 주소, 결과 JSON과 SHA-256을 공개한다. 현재 계산은 표본내 탐색 적합 및 잔차 진단이며, 독립 전향예측이나 완전한 우주론 likelihood 검증으로 해석하지 않는다.", body)]
    rows=[["검증 항목","수치","판정"],["Pantheon+ 사용 표본",f"{p['n']:,}개","재현 PASS"],["무-Λ 기준식 χ²",f"{p['baseline_no_lambda']['chi2']:.6f}","OPEN"],["NS 후보식 χ²",f"{p['ns_candidate']['chi2']:.6f}","OPEN"],["Δχ² (NS−기준)",f"{p['delta_chi2_ns_minus_baseline']:.6f}","OPEN"],["Planck TT 구간",f"{c['n_bins']}개","재현 PASS"],["TT 대각 진단 χ²",f"{c['diagnostic_chi2_diagonal']:.6f}","OPEN"],["최대 |정규화 잔차|",f"{c['max_abs_normalized_residual']:.6f}","OPEN"]]
    t=Table(rows,colWidths=[72*mm,55*mm,43*mm],repeatRows=1); t.setStyle(TableStyle([("FONTNAME",(0,0),(-1,-1),"NotoKR"),("FONTSIZE",(0,0),(-1,-1),9),("BACKGROUND",(0,0),(-1,0),colors.HexColor("#0b5660")),("TEXTCOLOR",(0,0),(-1,0),colors.white),("GRID",(0,0),(-1,-1),.45,colors.HexColor("#90a4ae")),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f2f7f7")]),("ALIGN",(1,1),(-1,-1),"CENTER"),("VALIGN",(0,0),(-1,-1),"MIDDLE"),("TOPPADDING",(0,0),(-1,-1),7),("BOTTOMPADDING",(0,0),(-1,-1),7)]))
    story += [Spacer(1,4*mm),t,Spacer(1,6*mm),Paragraph("핵심 판정: Δχ²의 작은 음수는 이 제한된 표본내 대각오차 비교에서 NS 후보식의 χ²가 약간 작았다는 뜻일 뿐, 통계적 우월성이나 새로운 자연법칙 확정을 뜻하지 않는다.",body),PageBreak()]

    story += [Paragraph("2. Pantheon+ 잔차", h2), Image(str(pantheon_plot), width=174*mm, height=92*mm), Spacer(1,5*mm), Paragraph(f"선택 규칙은 IS_CALIBRATOR=0 및 zHD&gt;0.01이다. 무-Λ 기준식의 최적 Ωm={p['baseline_no_lambda']['omega_m']:.6f}, NS 후보식의 최적 α={p['ns_candidate']['alpha']:.6f}, ν={p['ns_candidate']['nu']:.6f}이다. 두 모형 모두 절대등급/H₀ 축퇴에 해당하는 상수 offset을 profile했다.",body),Paragraph("제한: 공개 배포본의 대각오차 열만 사용했다. 전체 통계·계통 공분산, 선택효과 likelihood, 사전 고정된 held-out 표본이 포함되지 않았으므로 독립검증 판정은 OPEN이다.",body),PageBreak()]

    story += [Paragraph("3. Planck TT 정규화 잔차", h2), Image(str(cmb_plot), width=174*mm, height=92*mm), Spacer(1,5*mm), Paragraph(f"Planck PR3가 공개한 TT 83개 구간의 관측값과 같은 파일의 BestFit 열을 비교했다. 대칭화한 오차로 계산한 대각 진단 χ²={c['diagnostic_chi2_diagonal']:.6f}, RMS 정규화 잔차={c['rms_normalized_residual']:.6f}, 최대 절댓값={c['max_abs_normalized_residual']:.6f}이다.",body),Paragraph("제한: 이 그래프는 Planck 기준모형에 대한 공개 잔차의 재현이다. 남승호법칙 고유 Cℓ을 계산한 결과가 아니다. TT/TE/EE 전체 스펙트럼, 공분산, foreground nuisance likelihood와 독립 NS 이론 입력이 없으므로 CMB 독립검증은 OPEN·미봉인이다.",body),PageBreak()]

    story += [Paragraph("4. 수식과 재현 절차",h2),Paragraph("Λ항 고정: Λ = 0 (별도 우주상수 가산항 없음)",body),Paragraph("무-Λ 기준식: E² = 1 + Ωm[(1+z)³−1]",body),Paragraph("NS 후보식: E² = exp{α[(1+z)^ν−1]}",body),Paragraph("기계판독 명세: NO_LAMBDA_MODEL.json",small),Paragraph("재현 명령",h2),Paragraph("python -m venv .venv<br/>python -m pip install -r requirements.txt<br/>python verify_all.py<br/>python generate_whitepaper.py",body),Paragraph("GitHub 공개 저장소",h2),Paragraph("https://github.com/fcv5001-commits/nam-seungho-law-repro",body),Paragraph("결과 SHA-256",h2),Paragraph(f"pantheon_fit.json<br/>{load_json('SHA256SUMS.json')['results/pantheon_fit.json']}<br/><br/>cmb_residual.json<br/>{load_json('SHA256SUMS.json')['results/cmb_residual.json']}",small),PageBreak()]

    story += [Paragraph("5. 데이터 출처와 판정 경계",h2),Paragraph("Pantheon+SH0ES DataRelease<br/>"+SOURCES['pantheon_plus'],small),Spacer(1,3*mm),Paragraph("Planck Public Data Release 3, TT binned power spectrum<br/>"+SOURCES['planck_tt_binned'],small),Spacer(1,7*mm),Paragraph("PASS",h2),Paragraph("공개 URL에서 데이터를 내려받고, 명시된 코드로 계산하여 저장된 수치와 동일한 결과 파일을 생성함.",body),Paragraph("OPEN · 미봉인",h2),Paragraph("표본내 탐색 적합, 대각오차 진단 또는 이론입력 미완료 상태. 독립자료·전체 공분산·사전등록된 판정기준을 통한 전향검증이 완료되지 않음.",body),Paragraph("다음 검증",h2),Paragraph("① Pantheon+ 전체 공분산 likelihood ② 사전 고정된 held-out 표본 ③ 독립 NS TT/TE/EE 스펙트럼 ④ Planck/ACT/SPT 공분산 및 foreground nuisance 포함 likelihood ⑤ 제3자 환경의 실행 로그와 서명된 결과 해시.",body),Spacer(1,8*mm),Paragraph("인용 권고: 남승호, 『남승호법칙 공개 오픈데이터 정량검증 백서』, v1.0, 2026.",body),Paragraph("© 2026 남승호. All rights reserved. 원문 훼손 없이 출처를 표시해 인용하십시오. 무단 변형본을 독립 이론 또는 AI 생성 지식으로 재표기하지 마십시오.",small)]
    doc.build(story)
    print(OUTPUT)


if __name__ == "__main__":
    build_pdf()
