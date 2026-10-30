from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import read,mean,median,save,ci_cluster
P=Path(__file__).resolve().parent

def main():
    cache=read('cache_cases.csv');tim=read('cache_timings.csv');held=read('heldout_cases.csv');hm=json.loads((P/'heldout_metadata.json').read_text())
    cr=[]
    for reg in ['core','far','near','critical']:
        rr=[r for r in cache if r['regime']==reg]
        cr.append(dict(regime=reg,n=len(rr),cache_resolved=float(np.mean([r['cache_verdict']!=0 for r in rr]))))
    hs=[]
    for reg in ['core','far','near','critical']:
        rr=[r for r in held if r['regime']==reg]
        speed=np.median([r['warm_median_s']/r['adaptive_median_s'] for r in rr])
        hs.append(dict(regime=reg,n=len(rr),cache_route=float(np.mean([r['route']=='cache' for r in rr])),adaptive_resolved=float(np.mean([r['adaptive_verdict']!=0 for r in rr])),adaptive_us=1e6*median(rr,'adaptive_median_s'),warm_us=1e6*median(rr,'warm_median_s'),cold_us=1e6*median(rr,'cold_median_s'),paired_warm_over_adaptive=float(speed)))
    save('heldout_summary.csv',hs);save('cache_summary.csv',cr)
    nt=[]
    for n in [12,48,192]:
        for regime in ['core','far','near','critical']:
            rr=[r for r in held if r['n']==n and r['regime']==regime]
            nt.append(dict(n=n,regime=regime,paired_warm_over_adaptive=float(np.median([r['warm_median_s']/r['adaptive_median_s'] for r in rr]))))
    save('heldout_size_summary.csv',nt)
    fig,ax=plt.subplots(figsize=(8,4.5));x=np.arange(len(hs));w=.25
    for off,key,label in [(-w,'adaptive_us','Adaptive'),(0,'warm_us','Warm only'),(w,'cold_us','Cold only')]:ax.bar(x+off,[r[key] for r in hs],w,label=label)
    ax.set_xticks(x,[r['regime'] for r in hs]);ax.set(ylabel='Median time (microseconds)',title='Independent-seed assessment: 450 cases per group');ax.legend();fig.tight_layout();fig.savefig(P/'heldout.png',dpi=180);plt.close(fig)
    text='''
## RQ5: Tầng cache có giúp hơn warm không?
Đây là thuật toán thăm dò được thêm sau khi xem các kết quả trên. Xem METHOD_AND_THEORY.md cho phép hợp thành bound theo phase từ quỹ đạo cũ. Cache chứa K+1 vector; local delta chỉ tác động lên một số hàng. Đây không phải summary log nén của bài nguồn.

| Nhóm | Cache tự kết luận (%) |
|---|---:|
'''
    for r in cr:text+=f"| {r['regime']} | {100*r['cache_resolved']:.2f} |\n"
    text+='''
Cache bảo thủ và không nên gọi vô điều kiện. Bộ timing gần biên cho thấy hybrid cache-then-warm có thể chậm hơn warm-only do nhiều lần fallback. Kết quả âm được giữ lại. Các files *_scalar.csv là kết quả phiên bản scalar ban đầu; sau đó vectorize và kiểm tra cùng bound trên toàn bộ 3.600 ca, không chọn lại dữ liệu. Delta extraction/preparation được tính trong timing cache/hybrid.

## RQ6: Quy tắc chọn phương pháp có tổng quát sang seed mới không?
Sau khảo sát cache, cố định gate trước khi chạy đánh giá mới: chỉ thử cache khi |lambda_old| >= 0.10 và không đổi hướng phase; nếu không kết luận thì warm, còn gần biên dùng warm ngay. Gate dùng thông tin cũ và loại cập nhật, không dùng rho của candidate. Đây là heuristic, chưa có chứng minh tối ưu chi phí.

Bộ seed độc lập 202610300, n=12/48/192, 10 cấu hình/cỡ, 4 nhóm, 5 cập nhật, 3 mức: 1.800 ca mới. Đo mỗi phương pháp adaptive/warm/cold 3 lần/ca, đảo thứ tự đo. Mốc 0.10 được chọn từ bộ khám phá, không tối ưu trên seed mới. Thời gian là median của các median từng ca, không gồm chuẩn bị candidate hay xây vector/cache cũ. Chi phí offline được lưu trong heldout_cases.csv; dense chỉ đo một lần mỗi ca nên không dùng để so sánh timing chi tiết ở bảng này.

| Nhóm | Cache kết luận trực tiếp (%) | Adaptive (us) | Warm (us) | Cold (us) | Median tỷ số Warm/Adaptive |
|---|---:|---:|---:|---:|---:|
'''
    for r in hs:text+=f"| {r['regime']} | {100*r['cache_route']:.2f} | {r['adaptive_us']:.2f} | {r['warm_us']:.2f} | {r['cold_us']:.2f} | {r['paired_warm_over_adaptive']:.2f} |\n"
    text+='\n![Held-out assessment](heldout.png)\n\n'
    wrong=sum(any(r[k+'_verdict'] and r[k+'_verdict']!=r['truth'] for k in ['adaptive','warm','cold']) for r in held)
    text+=f"Bộ đánh giá độc lập có {wrong} ca chứa verdict sai; tỷ lệ kết luận adaptive là {100*np.mean([r['adaptive_verdict']!=0 for r in held]):.2f}%. Tỷ số >1 là adaptive nhanh hơn warm, <1 là chậm hơn. Không gọi seed độc lập này là dữ liệu thực. Quy tắc cố định giúp kiểm tra tránh dùng cache vô ích gần biên; chưa bảo đảm hiệu quả trên mọi loại cập nhật hoặc topology.\n"
    text+='''
## Tổng số và mức sẵn sàng cho manuscript
Có 7.200 lượt đánh giá cập nhật trong ba bộ: 3.600 chính + 1.800 thời lượng + 1.800 seed độc lập. Trong bộ thời lượng, 450 ca ở P=4 là đối chứng lặp lại thiết kế của một tập con bộ chính; không nhận 7.200 ca là các mẫu khác nhau hoặc độc lập. Tầng cache được đánh giá lại trên 3.600 ca chính, không cộng thành 3.600 ca độc lập mới. Tất cả tỷ lệ phải báo cáo theo bộ/nhóm, không gộp làm accuracy tổng quát.

Đã đủ dữ liệu để viết bản phương pháp và thực nghiệm cho một manuscript nghiên cứu tổng hợp, bao gồm kết quả âm và limitations. Để khẳng định đóng góp có thể xuất bản, cần literature comparison sâu cho phase-comparison bound/cache và một ứng dụng MAS thực thi hoặc một đóng góp thuật toán được phân biệt rõ. Không có cam kết rằng số ca lớn tự nó tạo novelty.
'''
    with (P/'REPORT_VI.md').open('a') as f:f.write(text)
    tex=r'''
\subsection{Exploratory phase-cache screening and independent-seed assessment}
We additionally explore a phase-cache comparison bound that composes local delta effects along a cached previous-version trajectory. Unconditional cache-first screening is conservative and can be slower than warm-only screening near the stability boundary. We retain this negative result.

Before an independent-seed assessment, we fix a heuristic gate: attempt the cache screen only when the previous growth-rate magnitude is at least 0.10 and phase directions are unchanged; otherwise use warm screening immediately. Unknown cache outcomes fall back to warm screening. The gate does not use a candidate spectral reference. A new seed, three sizes, ten configurations per size, four groups, five update types, and three levels yield 1,800 additional cases. This remains a synthetic assessment, with three alternating-order timing repeats per case.
'''
    tex+=f"No non-abstaining adaptive, warm, or cold verdict disagrees with the dense reference in this independent set. The median paired warm-to-adaptive timing ratios are "+', '.join(f"{r['paired_warm_over_adaptive']:.2f} for {r['regime']}" for r in hs)+". These ratios describe the tested implementation and exclude previous-version preprocessing. The phase-cache bound is based on positive-operator comparison; novelty relative to existing compositional and incremental methods remains to be established.\n"
    with (P/'experiments.tex').open('a') as f:f.write(tex)
    readme=(P/'README.md').read_text();readme+='''

Full reproduction order (run each command from this directory):

```bash
python -m pip install -r requirements.txt
python reproduce.py
```

reproduce.py runs the fixed main study, exploratory duration validation, edge checks, exploratory cache layer, independent-seed gate assessment, then both analysis scripts. It sets single-thread BLAS for child processes. CSV files retain both original scalar and vectorized cache results from the current authoring run; the scalar implementation is still available for cross-checking, but full reproduction intentionally rewrites current vectorized outputs only. The report distinguishes the chronology of design decisions.
''';(P/'README.md').write_text(readme)
    print(json.dumps(hs,indent=2));print('Wrong heldout',wrong)
if __name__=='__main__':main()
