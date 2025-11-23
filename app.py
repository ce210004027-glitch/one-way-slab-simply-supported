import streamlit as st
from io import BytesIO
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from one_way_slab_designer import run_design

st.set_page_config(page_title="210004027 - One way simply supported", layout='wide')
st.title("210004027 — One way simply supported")
st.markdown("**Manya Rajib Jain (210004027)**  
**Dr. Akshay Pratap Singh**  
**Civil Engineering, IIT Indore**")
st.write("---")

with st.form("inputs"):
    st.subheader("Geometry & Loads")
    col1, col2 = st.columns(2)
    with col1:
        Lx = st.number_input("Short span Lx (m)", value=3.0, format="%.3f")
        Ly = st.number_input("Long span Ly (m)", value=6.0, format="%.3f")
        h_mm = st.number_input("Overall slab thickness h (mm)", value=185.0, step=5.0)
        support_width_mm = st.number_input("Support width (mm)", value=230.0)
    with col2:
        live = st.number_input("Live load (kN/m^2)", value=4.0)
        floor_finish = st.number_input("Floor finish / external dead (kN/m^2)", value=1.5)
        fck = st.selectbox("Concrete grade fck (MPa)", options=[20,25,30], index=1)
        fy = st.selectbox("Steel grade fy (MPa)", options=[415,500], index=0)

    st.subheader("Reinforcement & Options")
    col3, col4 = st.columns(2)
    with col3:
        dia_bar = st.selectbox("Main bar diameter (mm)", options=[8,10,12,16], index=1)
        dist_bar = st.selectbox("Distribution bar diameter (mm)", options=[8], index=0)
    with col4:
        example_mode = st.checkbox("Example mode (tutorial defaults)", value=False)
        run_deflection = st.checkbox("Run IS deflection check (K-factor)", value=True)

    submitted = st.form_submit_button("Run Design")

if submitted:
    try:
        res = run_design(short_span_m=float(Lx), long_span_m=float(Ly), h_mm=float(h_mm),
                         fck=float(fck), fy=float(fy), external_dead_kNpm2=float(floor_finish),
                         live_kNpm2=float(live), dia_bar_mm=float(dia_bar), dist_bar_mm=float(dist_bar),
                         max_agg_mm=20.0)

        st.success("Design completed")
        st.subheader("Summary")
        st.write(f"Short span L = {res.short_span_m} m — Effective depth d = {res.d_mm:.1f} mm")
        st.write(f"Factored load w_u = {res.wu_kNpm:.3f} kN/m — Design moment M_u = {res.Mu_kNm_per_m:.3f} kN·m/m")
        st.write(f"Ast required = {res.Ast_required_mm2_per_m:.1f} mm^2/m — Ast provided = {res.Ast_provided_mm2_per_m:.1f} mm^2/m")
        st.write(f"Main bars: {int(res.dia_bar_mm)} mm @ {int(res.spacing_tension_mm)} mm c/c")
        st.write(f"Distribution: 8 mm @ {int(res.dist_spacing_mm)} mm c/c")
        st.write(f"Shear: V_u = {res.Vu_kN:.2f} kN, V_c = {res.Vc_kN:.2f} kN — shear_ok = {res.shear_ok}")
        st.write(f"L/d = {res.L_over_d:.2f}")

        if res.notes:
            st.subheader("Notes & warnings")
            for k, v in res.notes.items():
                st.warning(f"{k}: {v}")

        def generate_pdf(res):
            buffer = BytesIO()
            c = canvas.Canvas(buffer, pagesize=A4)
            text = c.beginText(40, 800)
            text.setFont("Helvetica", 10)
            text.textLine("One-way slab design report")
            text.textLine("")
            text.textLine(f"Short span L = {res.short_span_m} m")
            text.textLine(f"Overall depth h = {res.h_mm} mm, effective depth d = {res.d_mm:.1f} mm")
            text.textLine(f"M_u = {res.Mu_kNm_per_m:.3f} kN.m/m")
            text.textLine(f"Ast required = {res.Ast_required_mm2_per_m:.1f} mm2/m")
            text.textLine(f"Main: {int(res.dia_bar_mm)} mm @ {int(res.spacing_tension_mm)} mm c/c")
            text.textLine(f"Distribution: 8 mm @ {int(res.dist_spacing_mm)} mm c/c")
            c.drawText(text)
            c.showPage()
            c.save()
            buffer.seek(0)
            return buffer

        pdf_bytes = generate_pdf(res)
        st.download_button("Download PDF report", data=pdf_bytes, file_name="slab_design_report.pdf", mime='application/pdf')

    except Exception as e:
        st.error(f"Design failed: {e}")
