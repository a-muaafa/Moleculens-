import numpy as np
import pandas as pd
import plotly.graph_objects as go


def build_properties_chart(descriptors: dict) -> tuple[go.Figure, str]:
    labels = [
        "الوزن الجزيئي (من 500 Da)",
        "معامل التوزع LogP (من 5.0)",
        "مانحات الهيدروجين HBD (من 5)",
        "مستقبلات الهيدروجين HBA (من 10)",
        "الروابط الدوارة (من 10)",
        "المساحة القطبية TPSA (من 140 Å²)",
    ]
    raw_vals = [
        descriptors["molecular_weight"],
        descriptors["logp"],
        descriptors["hbd"],
        descriptors["hba"],
        descriptors["rotatable_bonds"],
        descriptors["tpsa"],
    ]
    thresholds = [500.0, 5.0, 5.0, 10.0, 10.0, 140.0]
    pct_vals = [max(0.0, (v / t) * 100.0) for v, t in zip(raw_vals, thresholds)]
    colors = ["#15803d" if p <= 100.0 else "#b91c1c" for p in pct_vals]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=labels,
            x=pct_vals,
            orientation="h",
            marker_color=colors,
            text=[f"{v} ({p:.0f}%)" for v, p in zip(raw_vals, pct_vals)],
            textposition="auto",
            name="نسبة الخاصية إلى الحد المرجعي",
        )
    )
    fig.add_vline(
        x=100.0,
        line_dash="dash",
        line_color="#0f3460",
        annotation_text="الحد المرجعي الاسترشادي (100%)",
        annotation_position="top right",
    )
    fig.update_layout(
        title="مقارنة الخصائص الفيزيائية الكيميائية للمركب بالحدود المرجعية الاسترشادية (Lipinski / Veber)",
        xaxis_title="النسبة المئوية من الحد المرجعي (%)",
        yaxis_title="الخاصية الكيميائية",
        height=320,
        margin=dict(l=20, r=20, t=45, b=20),
    )
    explanation = (
        "تفسير الرسم: يوضح الشريط الأفقي موقع كل خاصية كيميائية نسبة إلى الحد الأعلى الاسترشادي للامتصاص الفموي (100%). "
        "الأشرطة الخضراء تعني بقاء الخاصية ضمن المدى المرجعي، بينما تشير الأشرطة الحمراء إلى تجاوز الحد الاسترشادي."
    )
    return fig, explanation


def build_published_activity_distribution_chart(
    clean_df: pd.DataFrame,
    highlight_pic50: float | None = None,
    highlight_label: str = "قيمة المركب",
) -> tuple[go.Figure, str]:
    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=clean_df["pIC50"],
            nbinsx=28,
            marker_color="#0f3460",
            opacity=0.82,
            name="قياسات ChEMBL المنشورة (pIC50)",
        )
    )
    fig.add_vline(
        x=6.0,
        line_dash="dash",
        line_color="#b45309",
        annotation_text="عتبة النشاط الاسترشادية (pIC50 = 6.0 / IC50 = 1000 nM)",
        annotation_position="top left",
    )
    if highlight_pic50 is not None:
        fig.add_vline(
            x=float(highlight_pic50),
            line_width=3,
            line_color="#15803d",
            annotation_text=f"{highlight_label} ({highlight_pic50:.2f})",
            annotation_position="top right",
        )
    fig.update_layout(
        title="توزيع القياسات المخبرية المنشورة في ChEMBL لهدف DPP4 البشري (على مقياس pIC50)",
        xaxis_title="قيمة pIC50 المخبرية (كلما ارتفعت القيمة زادت قوة التثبيط وانخفض تركيز IC50)",
        yaxis_title="عدد المركبات الفريدة",
        height=310,
        margin=dict(l=20, r=20, t=45, b=20),
    )
    explanation = (
        "تفسير الرسم: يعرض المدرج التكراري توزيع النشاط الحيوي الحقيقي المقاس مخبريا لمركبات ChEMBL ضد إنزيم DPP4 البشري. "
        "تمثل القيم على يمين الخط المتقطع (pIC50 >= 6.0) المركبات ذات تركيز تثبيط نصفي IC50 يبلغ 1000 nM أو أقل."
    )
    return fig, explanation


def build_pic50_gauge(pic50_val: float) -> tuple[go.Figure, str]:
    clamped_val = float(np.clip(pic50_val, 3.0, 10.0))
    fig = go.Figure(
        go.Indicator(
            mode="gauge+number",
            value=clamped_val,
            number={"suffix": " pIC50", "valueformat": ".2f"},
            title={"text": "مقياس التنبؤ الحاسوبي اللوغاريتمي (pIC50 من 3.0 إلى 10.0)"},
            gauge={
                "axis": {"range": [3.0, 10.0], "tickvals": [3, 4, 5, 6, 7, 8, 9, 10]},
                "bar": {"color": "#0f3460"},
                "steps": [
                    {"range": [3.0, 6.0], "color": "#fecaca"},
                    {"range": [6.0, 7.0], "color": "#fde68a"},
                    {"range": [7.0, 10.0], "color": "#bbf7d0"},
                ],
                "threshold": {
                    "line": {"color": "#0b192c", "width": 4},
                    "thickness": 0.8,
                    "value": clamped_val,
                },
            },
        )
    )
    fig.update_layout(height=275, margin=dict(l=25, r=25, t=55, b=15))
    explanation = (
        "تفسير المقياس: يمثل pIC50 اللوغاريتم السالب لتركيز IC50 بالمولار؛ المنطقة الحمراء (3-6) تقابل IC50 أعلى من 1000 nM، "
        "والمنطقة الصفراء (6-7) تقابل 100 إلى 1000 nM، والمنطقة الخضراء (7-10) تقابل أقل من 100 nM."
    )
    return fig, explanation


def build_similarity_chart(all_sims: np.ndarray, query_max_sim: float) -> tuple[go.Figure, str]:
    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=all_sims,
            nbinsx=30,
            marker_color="#0d5c63",
            opacity=0.82,
            name="مركبات التدريب",
        )
    )
    fig.add_vline(
        x=0.30,
        line_dash="dash",
        line_color="#b91c1c",
        annotation_text="حد خارج النطاق (0.30)",
        annotation_position="top left",
    )
    fig.add_vline(
        x=0.60,
        line_dash="dash",
        line_color="#15803d",
        annotation_text="حد داخل النطاق (0.60)",
        annotation_position="top left",
    )
    fig.add_vline(
        x=query_max_sim,
        line_width=3,
        line_color="#0f3460",
        annotation_text=f"أقصى تشابه للمركب ({query_max_sim:.2f})",
        annotation_position="top right",
    )
    fig.update_layout(
        title="توزيع تشابه Tanimoto للمركب المدخل مقارنة ببصمات مركبات التدريب فقط",
        xaxis_title="معامل تشابه Tanimoto مع مركبات التدريب",
        yaxis_title="عدد مركبات التدريب",
        height=310,
        margin=dict(l=20, r=20, t=45, b=20),
    )
    explanation = (
        "تفسير الرسم: يقارن الرسم بصمة المركب المدخل بجميع مركبات التدريب لتحديد موقعه بالنسبة لنطاق تطبيق النموذج: "
        "أقل من 0.30 (خارج النطاق)، بين 0.30 و 0.59 (حدي)، و 0.60 فأعلى (داخل النطاق)."
    )
    return fig, explanation


def build_nearest_compounds_chart(nearest_df: pd.DataFrame) -> tuple[go.Figure, str]:
    ids = nearest_df["molecule_chembl_id"].astype(str).tolist()
    sims = nearest_df["tanimoto"].astype(float).tolist()
    pic50s = nearest_df["pIC50"].astype(float).tolist()

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=ids,
            y=sims,
            name="معامل التشابه (Tanimoto)",
            marker_color="#0d5c63",
            text=[f"{s:.2f}" for s in sims],
            textposition="auto",
            yaxis="y1",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=ids,
            y=pic50s,
            name="النشاط المخبري (pIC50)",
            mode="lines+markers",
            marker=dict(size=10, color="#0f3460"),
            yaxis="y2",
        )
    )
    fig.update_layout(
        title="أقرب 5 مركبات في بيانات التدريب: درجة التشابه الهيكلي مقابل النشاط المخبري المنشور",
        xaxis_title="معرف المركب في بيانات التدريب (ChEMBL ID)",
        yaxis=dict(title="معامل التشابه (Tanimoto)", range=[0, 1.05]),
        yaxis2=dict(title="النشاط المخبري (pIC50)", overlaying="y", side="left", range=[3, 10]),
        height=320,
        margin=dict(l=40, r=40, t=45, b=20),
    )
    explanation = (
        "تفسير الرسم: تعرض الأعمدة درجة تشابه Tanimoto لأقرب 5 مركبات في مجموعة التدريب، بينما يوضح الخط النقطي "
        "قيمة النشاط المخبري الحقيقي (pIC50) المسجلة لكل منها في ChEMBL."
    )
    return fig, explanation


def build_split_comparison_chart(metrics: dict) -> tuple[go.Figure, str]:
    rand = metrics["random_split_comparison"]
    metric_names = ["MAE (أقل أفضل)", "RMSE (أقل أفضل)", "R² (أعلى أفضل)", "ROC-AUC (للتصنيف)"]
    scaf_vals = [metrics["mae_pic50"], metrics["rmse_pic50"], metrics["r2_pic50"], metrics["roc_auc"]]
    rand_vals = [rand["mae_pic50"], rand["rmse_pic50"], rand["r2_pic50"], rand["roc_auc"]]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=metric_names,
            y=scaf_vals,
            name="Scaffold Split (اختبار هياكل جديدة)",
            marker_color="#0f3460",
            text=[f"{v:.3f}" for v in scaf_vals],
            textposition="auto",
        )
    )
    fig.add_trace(
        go.Bar(
            x=metric_names,
            y=rand_vals,
            name="Random Split (تقسيم عشوائي)",
            marker_color="#b45309",
            text=[f"{v:.3f}" for v in rand_vals],
            textposition="auto",
        )
    )
    fig.update_layout(
        barmode="group",
        title="مقارنة أداء النموذج بين التقسيم الهيكلي (Scaffold Split) والتقسيم العشوائي (Random Split)",
        xaxis_title="مقياس التقييم",
        yaxis_title="القيمة",
        height=340,
        margin=dict(l=20, r=20, t=45, b=20),
    )
    explanation = (
        "تفسير الرسم: يظهر التقسيم العشوائي أرقاما أكثر تفاؤلا (خطأ أقل و R² أعلى) بسبب وجود هياكل أساسية مشتركة بين التدريب والاختبار، "
        "بينما يعكس Scaffold Split الأداء الواقعي للنموذج عند مواجهة عائلات كيميائية جديدة تماما."
    )
    return fig, explanation


def build_test_scatter_chart(test_df: pd.DataFrame) -> tuple[go.Figure, str]:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=test_df["true_pIC50"],
            y=test_df["pred_pIC50"],
            mode="markers",
            marker=dict(
                size=9,
                color=test_df["max_train_tanimoto"],
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title="Tanimoto"),
            ),
            text=(
                test_df["molecule_chembl_id"].astype(str)
                + " | خطأ مطلق: "
                + test_df["abs_error_pIC50"].round(2).astype(str)
                + " | Tanimoto: "
                + test_df["max_train_tanimoto"].round(2).astype(str)
            ),
            hoverinfo="text",
            name="مركبات الاختبار",
        )
    )
    min_v = float(min(test_df["true_pIC50"].min(), test_df["pred_pIC50"].min()) - 0.3)
    max_v = float(max(test_df["true_pIC50"].max(), test_df["pred_pIC50"].max()) + 0.3)
    fig.add_trace(
        go.Scatter(
            x=[min_v, max_v],
            y=[min_v, max_v],
            mode="lines",
            line=dict(color="#b91c1c", dash="dash"),
            name="خط التطابق المثالي",
        )
    )
    fig.update_layout(
        title="القيم المخبرية الحقيقية مقابل المتوقعة حاسوبيا على مجموعة اختبار Scaffold Split",
        xaxis_title="قيمة pIC50 المخبرية الحقيقية (ChEMBL)",
        yaxis_title="قيمة pIC50 المتوقعة حاسوبيا (Random Forest)",
        height=370,
        margin=dict(l=20, r=20, t=45, b=20),
    )
    explanation = (
        "تفسير الرسم: تمثل كل نقطة مركبا في مجموعة الاختبار ذات الهياكل غير المرئية في التدريب، ملونة بدرجة تشابهها مع التدريب. "
        "النقاط القريبة من الخط المتقطع تمثل تنبؤات دقيقة، بينما تكشف النقاط البعيدة عن حدود التعميم وقفزات النشاط."
    )
    return fig, explanation

