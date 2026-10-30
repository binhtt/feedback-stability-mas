"""Generate manuscript-facing tables, figures and a candid Vietnamese report."""
from pathlib import Path
import csv,json,platform,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent

def read(name):
    rows=list(csv.DictReader((P/name).open()))
    for r in rows:
        for k,v in r.items():
            try:r[k]=float(v)
            except ValueError:pass
    return rows

def mean(rs,key):return float(np.mean([r[key] for r in rs]))
def median(rs,key):return float(np.median([r[key] for r in rs]))
def save(name,rows):
    with (P/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)

def ci_cluster(rs,key,B=2000):
    """Stratified paired cluster bootstrap: config is a cluster, strata are n."""
    rng=np.random.default_rng(202610291);sizes=sorted({r['n'] for r in rs});strata=[]
    for n in sizes:
        rr=[r for r in rs if r['n']==n];reps=sorted({r['rep'] for r in rr})
        vals=[]
        for rep in reps:
            cr=[r for r in rr if r['rep']==rep]
            vals.append(float(np.mean([key(r) if callable(key) else r[key] for r in cr])))
        strata.append(np.array(vals))
    draws=np.zeros(B)
    for vals in strata:draws+=vals[rng.integers(0,len(vals),size=(B,len(vals)))].mean(axis=1)/len(strata)
    return np.quantile(draws,[.025,.975]).tolist()

def main():
    rows=read('cases.csv');ts=read('timings.csv');ext=read('duration_sensitivity.csv');base=read('baselines.csv');m=json.loads((P/'metadata.json').read_text());em=json.loads((P/'extensions_metadata.json').read_text())
    stats=[]
    for reg in ['core','far','near','critical']:
        rr=[r for r in rows if r['regime']==reg]
        for budget in [1,2,4,8,16]:
            warm=mean(rr,f'warm_resolved_{budget}');cold=mean(rr,f'cold_resolved_{budget}')
            ci=ci_cluster(rr,lambda r:r[f'warm_resolved_{budget}']-r[f'cold_resolved_{budget}'])
            stats.append(dict(regime=reg,budget=budget,n=len(rr),warm_rate=warm,cold_rate=cold,difference=warm-cold,paired_cluster_ci_low=ci[0],paired_cluster_ci_high=ci[1]))
    save('resolution_summary.csv',stats)
    timing=[]
    for n in m['sizes']:
        rr=[r for r in ts if r['n']==n]
        timing.append(dict(n=n,n_scenarios=len(rr),dense_ms=1000*median(rr,'dense_median_s'),warm_ms=1000*median(rr,'warm_median_s'),cold_ms=1000*median(rr,'cold_median_s'),arpack_ms=1000*median(rr,'arpack_median_s'),median_paired_dense_over_warm=float(np.median([r['dense_median_s']/r['warm_median_s'] for r in rr])),median_paired_cold_over_warm=float(np.median([r['cold_median_s']/r['warm_median_s'] for r in rr])),max_arpack_error=max(r['arpack_max_rho_error'] for r in rr)))
    save('timing_summary.csv',timing)
    durations=[]
    for h in [.005,.025,.125,.5]:
        rr=[r for r in ext if r['h']==h]
        ci=ci_cluster(rr,lambda r:r['warm_resolved_16']-r['cold_resolved_16'])
        durations.append(dict(h=h,period=32*h,n=len(rr),warm1=mean(rr,'warm_resolved_1'),cold1=mean(rr,'cold_resolved_1'),warm16=mean(rr,'warm_resolved_16'),cold16=mean(rr,'cold_resolved_16'),warm_steps=mean(rr,'warm_steps'),cold_steps=mean(rr,'cold_steps'),difference_ci_low=ci[0],difference_ci_high=ci[1]))
    save('duration_summary.csv',durations)
    # Failure and bound checks are assertions in study and extensions, retained logs too.
    wrong=[r for r in rows if any(r[x+'_verdict'] and r[x+'_verdict']!=r['truth'] for x in ['warm','cold'])]
    assert not wrong
    # Main resolution curves.
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),sharey=True)
    for ax,reg in zip(axs,['core','near','critical']):
        rr=[r for r in stats if r['regime']==reg]
        ax.plot([r['budget'] for r in rr],[100*r['warm_rate'] for r in rr],'-o',label='Previous-version vector')
        ax.plot([r['budget'] for r in rr],[100*r['cold_rate'] for r in rr],'--s',label='Ones vector')
        ax.set(xlabel='Period passes (budget)',title=reg,ylim=(-3,103));ax.set_xscale('log',base=2);ax.set_xticks([1,2,4,8,16],['1','2','4','8','16']);ax.grid(alpha=.2)
    axs[0].set_ylabel('Resolved decisions (%)');axs[-1].legend(fontsize=8);fig.tight_layout();fig.savefig(P/'resolution.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.5))
    for key,label in [('dense_ms','Dense product + eig'),('warm_ms','Warm CW screen'),('cold_ms','Cold CW screen'),('arpack_ms','Matrix-free ARPACK (point estimate)')]:
        ax.loglog([r['n'] for r in timing],[r[key] for r in timing],'-o',label=label)
    ax.set(xlabel='Number of agents',ylabel='Median online time (ms)',title='Timing: 9 cases/size, 7 repeats/case, one BLAS thread');ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(P/'timing.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.5));x=np.arange(4);width=.35
    ax.bar(x-width/2,[100*r['warm16'] for r in durations],width,label='Previous-version vector');ax.bar(x+width/2,[100*r['cold16'] for r in durations],width,label='Ones vector');ax.set_xticks(x,[str(r['period']) for r in durations]);ax.set(xlabel='Period duration P (exploratory extension)',ylabel='Resolved within 16 passes (%)',ylim=(0,108));ax.legend();fig.tight_layout();fig.savefig(P/'duration.png',dpi=180);plt.close(fig)
    text='''# Expanded experimental study

## Kết luận chính
Đã chạy 3.600 ca chính và 1.800 ca khảo sát thời lượng tương tác; mã nguồn và dữ liệu đầy đủ đi kèm. Tái sử dụng vector dương từ phiên bản cũ cải thiện tỷ lệ kết luận sớm. Trong bộ chính, cả warm và cold đều kết luận hết trước hoặc tại 16 lượt, không có kết luận sai so với dense reference. Trong bộ mở rộng chu kỳ ngắn, có ca chưa kết luận; unknown không được coi là stable.

Đây là nghiên cứu tổng hợp trên mô hình tuyến tính dương, không phải bằng chứng thực nghiệm về LLM hoặc MAS production. Thuật toán dùng bất đẳng thức Collatz–Wielandt và power iteration đã biết; không nhận các công cụ toán này là đóng góp mới. Kết quả chưa đủ để xác nhận novelty của toàn bài.

## Thiết kế cố định của bộ chính
- Seed 202610290, khác pilot; không lọc seed hoặc thay thế cấu hình theo kết quả.
- n = 6, 12, 24, 48, 96, 192; 10 cấu hình mỗi cỡ: tổng 60 cấu hình độc lập theo generator.
- Hai nhóm cùng cỡ. F,G có entries Uniform(0.1,1), chuẩn hóa theo trung bình tổng hàng; decay Uniform(0.8,1.2).
- 16 pulse cân bằng X/Y và 16 khoảng nghỉ, h=0.125, P=4.
- core: beta Uniform(1,6), không hiệu chỉnh trạng thái ổn định.
- far/near/critical: hiệu chỉnh beta để lambda gốc lần lượt -0.20/-0.02/-0.002. Đây là stress test có chủ đích, báo cáo riêng core.
- 5 loại cập nhật x 3 mức x 4 nhóm x 60 cấu hình = 3.600 ca. Các ca trên cùng cấu hình có phụ thuộc thống kê.
- weight: tăng hàng 0 F. decay: giảm decay agent 0. joint: tăng hàng 0 của F và G. redistribute: trộn hàng 0 của F với phiên bản dịch vòng, bảo toàn tổng hàng. order: 1/2/4 swaps Y sớm nhất với X muộn nhất, giữ inventory/duration/ma trận. Nhãn 5/10/20% chỉ có nghĩa phần trăm cho cập nhật tham số; với order là mức thao tác danh nghĩa.

## Phương pháp và baseline
Warm: lấy vector Perron dương của phiên bản cũ, truyền qua từng phase của phiên bản mới; dựng min/max tỷ số thành phần. Cold: cùng thuật toán nhưng vector đầu toàn 1. Tối đa 16 lượt, dừng ngay khi kết luận; guard float64 1e-10 quanh rho=1.

Dense reference: cùng công thức phase giải tích nhưng dựng đầy đủ Phi và tính eigenvalues. ARPACK: LinearOperator matrix-free, v0 toàn 1, tol=1e-10; là point estimator, không có cùng ngữ nghĩa abstention/bound với CW. Không gộp ARPACK thành phương pháp chứng nhận.

Preparation được ghi riêng. Online timing loại trừ lấy vector Perron của phiên bản cũ; chi phí offline đó được lưu trong baselines.csv. Chi phí warm không phải chi phí cold start. Độ phức tạp sparse/structured propagation: O(b K n²); dense reference O(K n³+n³). Cả warm và cold đều hưởng lợi từ matrix-free structure. Bản triển khai hiện tại vẫn chuẩn bị lại phase kernels, chưa cập nhật cache theo sparse delta; do đó chưa gọi là thuật toán incremental đầy đủ.

## RQ1: Tái sử dụng phiên bản cũ giúp kết luận sớm đến đâu?

Bảng tỷ lệ kết luận sau một lượt, mỗi nhóm 900 ca. Khoảng tin cậy cho chênh lệch dùng paired cluster bootstrap, phân tầng theo cỡ hệ, 2.000 lần, seed cố định; cluster là cấu hình. Đây là bất định trong generator tổng hợp, không phải tổng thể MAS thực tế.

| Nhóm | Warm (%) | Cold (%) | Chênh lệch (điểm %) | Bootstrap 95% CI |
|---|---:|---:|---:|---|
'''
    for reg in ['core','far','near','critical']:
        r=next(r for r in stats if r['regime']==reg and r['budget']==1)
        text+=f"| {reg} | {100*r['warm_rate']:.2f} | {100*r['cold_rate']:.2f} | {100*r['difference']:.2f} | [{100*r['paired_cluster_ci_low']:.2f}, {100*r['paired_cluster_ci_high']:.2f}] |\n"
    text+='\n![Resolution](resolution.png)\n\nỞ budget 4, cả hai giải quyết gần như toàn bộ bộ chính; ở budget 16 cả hai đạt 100%. Vì vậy, không dùng tốc độ so với dense để khẳng định riêng warm-start là nguồn của toàn bộ cải thiện.\n'
    text+='''
## RQ2: Chi phí theo cỡ hệ
9 kịch bản/cỡ hệ: 3 cấu hình đầu, near, mức 10%, các loại weight/joint/order; 7 lần đo mỗi phương pháp, đảo thứ tự đo luân phiên. Các số ms là median qua 9 median từng ca. Bảng tốc độ trong CSV dùng median của các tỷ số ghép cặp, khác với tỷ số của median. Thời gian chỉ có giá trị trên môi trường hiện tại; single-thread BLAS, không suy ra hiệu năng production.

| n | Dense (ms) | Warm (ms) | Cold (ms) | ARPACK (ms) |
|---|---:|---:|---:|---:|
'''
    for r in timing:text+=f"| {r['n']} | {r['dense_ms']:.4f} | {r['warm_ms']:.4f} | {r['cold_ms']:.4f} | {r['arpack_ms']:.4f} |\n"
    text+='\n![Timing](timing.png)\n\nCold có thể nhanh tương đương hoặc hơn warm ở một số cỡ lớn. Vì cập nhật một hàng chiếm tỷ lệ giảm theo n, cỡ hệ và độ mạnh toàn cục của cập nhật bị confound; chưa thể kết luận hiệu quả warm tăng đơn điệu theo n. Cần thêm cập nhật theo tỷ lệ agent cố định khi nghiên cứu riêng tác động của n.\n'
    text+='''
## RQ3: Chu kỳ ngắn có làm lợi ích warm rõ hơn không?
Đây là extension thăm dò được thiết kế sau khi xem bộ chính, không thay thế kết quả chính. Dùng n=24/96, 5 cấu hình đầu/cỡ, 3 stress margins, 5 cập nhật, 3 mức, 4 thời lượng: tổng 1.800 ca. Mỗi thời lượng có 450 ca. Beta được hiệu chỉnh lại bằng dense root để giữ cùng lambda gốc, nhưng thời lượng vẫn làm thay đổi kernel và khoảng cách phổ; không coi đây là thay đổi thuần một yếu tố vật lý.

| P | Warm kết luận/16 lượt (%) | Cold kết luận/16 lượt (%) | Số lượt TB warm | Số lượt TB cold |
|---|---:|---:|---:|---:|
'''
    for r in durations:text+=f"| {r['period']:.2f} | {100*r['warm16']:.2f} | {100*r['cold16']:.2f} | {r['warm_steps']:.2f} | {r['cold_steps']:.2f} |\n"
    text+='\n![Duration](duration.png)\n\nSố lượt trung bình bị chặn tại 16; đây không phải số lượt hội tụ thực của các ca unknown. Kết quả hỗ trợ dùng warm vector khi ngân sách thấp và mixing chậm, nhưng chưa phải chứng minh hội tụ nhanh hơn trên mọi hệ.\n'
    text+=f'''
## RQ4: Kiểm tra độ đúng
- 3.600 ca chính: không có verdict sai so với dense reference; bounds đối chiếu có tolerance 1e-8.
- 1.800 ca duration: không có verdict sai khi phương pháp đưa kết luận; unknown được giữ nguyên.
- {em['n_independent_checks']} kiểm tra phase giải tích và structured action so với dense expm triển khai độc lập.
- {em['n_interval_checks']} kiểm tra lại bằng mpmath interval arithmetic 60 chữ số trên n=6/12, critical margin, 5 mutation types, mức 10%. Mỗi candidate dùng vector dương sau 16 bước, rồi độc lập truyền interval để xác nhận dấu. Đây là interval verification trên tập con, không làm 5.400 verdict float64 thành certificate.
- 5 edge checks: stable/unstable/reducible/exact boundary, scaling của vector, monotonicity coupling.

Mô hình trong interval checks dùng chính các binary64 inputs như số thực chính xác; không xét sai số đo của dữ liệu. Bản thân vector khởi đầu có thể gần đúng, vì bất đẳng thức đúng cho mọi vector dương.

## Điều đã có và phần còn thiếu trước khi submit
Đã có: model rõ, lý thuyết sàng lọc, thuật toán thực thi, strong ablation warm-vs-cold, dense reference, ARPACK, 5.400 ca, phân tầng stress/core, repeated timing, uncertainty reporting, interval subset, artifact tái lập và văn bản phương pháp/thực nghiệm.

Chưa có: đóng góp mới được xác nhận so với literature incremental positive/switched systems; MAS thực thi có ánh xạ được kiểm tra; dữ liệu thực; history-dependent hoặc nonperiodic guarantees; thí nghiệm cập nhật tuần tự dài hạn; chứng nhận interval toàn bộ; kernel sparse/cache cập nhật theo delta; sửa lỗi hoặc tự động repair.

Nếu bài nhắm formal/algorithmic venue, cần đóng góp mới vượt warm-start CW cổ điển. Nếu nhắm applied MAS/software engineering, cần thêm case study MAS thực thi và chứng minh model phản ánh propagation của hệ đó. Không trình bày kết quả tổng hợp như đánh giá an toàn LLM. Các tỷ lệ là mô tả trong thiết kế generator, không tỷ lệ lỗi thực tế.

## Tái lập
1. Cài requirements.txt trong môi trường Python riêng.
2. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python study.py
3. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extensions.py
4. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python test_checks.py
5. python analyze.py

Thông tin môi trường, thiết kế và runtime có trong metadata.json/extensions_metadata.json. Không cần mạng, GPU hay LLM API sau cài dependency. Các bước ghi đè kết quả có chủ đích; lưu bản cũ nếu muốn so sánh. configs.json giữ đầy đủ ma trận và lịch sinh. Ảnh dùng tiếng Anh để có thể đưa vào manuscript.
'''
    (P/'REPORT_VI.md').write_text(text)
    # Ready-to-edit experiment prose; no inflated novelty, no fabricated significance.
    core=next(r for r in stats if r['regime']=='core' and r['budget']==1)
    near=next(r for r in stats if r['regime']=='near' and r['budget']==1)
    critical=next(r for r in stats if r['regime']=='critical' and r['budget']==1)
    d=durations[0]
    tex=r'''\section{Experimental Evaluation}
\label{sec:evaluation}
We evaluate whether reusing a positive vector from a previous system version improves the early resolution of feedback-stability queries. The experiments concern synthetic periodic positive linear systems; they do not evaluate the semantic correctness or safety of LLM agents.

\subsection{Design and baselines}
We generate ten configurations at each of six sizes, $n\in\{6,12,24,48,96,192\}$. Each configuration contains two equally sized populations, positive cross-population kernels, heterogeneous diagonal losses, and a balanced sixteen-pulse schedule interleaved with idle intervals. The period is $P=4$. We retain an uncalibrated core group with interaction intensity sampled uniformly from $[1,6]$, and separate stress groups whose original growth rates are calibrated to $-0.20$, $-0.02$, and $-0.002$. Five update types at three levels yield 3,600 cases. Updates sharing a configuration are paired observations rather than independent samples.

We compare a previous-version Perron-vector initialization (warm) with an all-ones initialization (cold), using the same Collatz--Wielandt screening procedure and a budget of sixteen period passes. Dense period-product construction followed by eigendecomposition provides a numerical reference. Matrix-free ARPACK provides an additional point-estimation timing baseline, not a bound-based verification procedure. Candidate preparation and previous-version preprocessing are recorded separately.

\subsection{Early resolution}
'''
    tex+=f"After one pass, warm and cold resolve {100*core['warm_rate']:.2f}\\% and {100*core['cold_rate']:.2f}\\% of core cases, respectively. The corresponding rates are {100*near['warm_rate']:.2f}\\% versus {100*near['cold_rate']:.2f}\\% in the near group and {100*critical['warm_rate']:.2f}\\% versus {100*critical['cold_rate']:.2f}\\% in the critical group. Both methods resolve all main-study cases within sixteen passes. No non-abstaining verdict disagrees with the dense numerical reference. Thus, the main benefit observed here is earlier resolution rather than a difference in eventual coverage.\n\n"
    tex+=r'''\subsection{Scalability and comparison with existing approaches}
Timing uses nine scenarios per size and seven repetitions per scenario, with alternating method order and single-threaded BLAS. Both warm and cold matrix-free screening avoid dense period-product construction. Cold screening can match or outperform warm screening on individual timing subsets; the improvement over dense recomputation must therefore not be attributed entirely to previous-version reuse. Previous-version eigenvector preprocessing is excluded from online timings and is separately reported. Local updates affect a decreasing fraction of agents as size increases, which limits causal interpretation of scaling trends.

\subsection{Exploratory duration sensitivity}
After inspecting the main study, we fix an additional duration-sensitivity experiment with 1,800 cases, retaining sizes 24 and 96 and five configurations per size. This extension is exploratory. Original growth-rate margins are recalibrated at each duration.
'''
    tex+=f"At period $P={d['period']:.2f}$, warm resolves {100*d['warm16']:.2f}\\% of cases within sixteen passes, compared with {100*d['cold16']:.2f}\\% for cold. Remaining cases are reported as unknown. These results indicate a benefit under a limited pass budget in the tested short-period regime, rather than universal superiority.\n\n"
    tex+=r'''\subsection{Validation and limitations}
We independently compare analytic phase propagation with dense matrix exponentials, test exact-boundary and reducible examples, and recheck thirty small-system candidates using 60-digit interval arithmetic. These interval checks validate a subset of supplied numerical instances; they do not certify every floating-point result. Paired cluster bootstrap intervals use configuration-level clusters stratified by size. Synthetic generation, calibrated stress margins, fixed periodic schedules, dense positive kernels, and the absence of an executable application case study limit external validity. The screening inequalities and power iteration are established mathematical tools; novelty must be assessed at the level of the evolution-analysis method and its application.
'''
    (P/'experiments.tex').write_text(tex)
    (P/'requirements.txt').write_text(f"numpy=={m['numpy']}\nscipy=={m['scipy']}\nmatplotlib=={matplotlib.__version__}\nmpmath=={em['mpmath']}\n")
    (P/'README.md').write_text('''# Feedback Stability of Evolving Multi-Agent Systems

Start with REPORT_VI.md (results and limitations), METHOD_AND_THEORY.md (model and exact-arithmetic bounds), and experiments.tex (editable manuscript text). Run study.py, extensions.py, test_checks.py, then analyze.py, in that order. Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before importing NumPy. Install requirements.txt first. The full main and exploratory designs and outputs are retained, including unknown cases. This is synthetic experimental evidence, not a verified LLM MAS or a novelty claim.
''')
    print('Main cases',len(rows),'extension',len(ext),'wrong main',len(wrong));print('ONE PASS',json.dumps([r for r in stats if r['budget']==1],indent=2));print('TIMING',json.dumps(timing,indent=2));print('DURATION',json.dumps(durations,indent=2))
if __name__=='__main__':main()
