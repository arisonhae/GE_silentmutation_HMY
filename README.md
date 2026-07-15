# GE_silentmutation_HMY

WT/Edited 시퀀스로부터 pegRNA를 설계하고, 각 RTT 안에서 MMR을 회피하는 **silent
co-edit**(동의 코돈 치환)을 찾아 prime editing 효율을 높이는 웹툴. reading frame과
exon 구간을 직접 지정할 수 있어, window에 인트론이 섞여 있어도(edit이 인트론에 있거나
exon 경계에 걸친 경우 포함) 인트론을 코돈으로 잘못 읽지 않고 exon 안에서만 silent
마커를 찾는다. Substitution과 1–3bp indel을 지원하며, 효율 예측은 DeepPrime 기반.

> **⚠️ 서열에 인트론이 섞여 있으면 exon-aware 체크박스를 켜고 exon 구간을 마킹할 것
> (꺼짐이 기본값).** 끄면 툴이 창 전체를 CDS로 가정해 인트론을 코돈으로 읽고, 인트론
> 염기에 "silent" 치환을 제안하게 된다. silent를 확인하는 window는 변이 기준 ±60bp라,
> 변이가 exon 가장자리 근처면 인트론이 딸려 들어오게 된다.

---

## 두 가지 사용법

### 1) 설치 없이 바로 사용

`silent_mutation_standalone.html` **파일을 브라우저에서 열기**
파이썬·conda·서버 불필요

- 첫 실행 때 파이썬 런타임(Pyodide)을 CDN에서 받으므로 **인터넷이 필요**함
  (10–20초, 이후엔 브라우저 캐시라 빠르게 실행 가능).
- 해당 파일 안에 `silent_mutation/` 코어가 동일하게 들어 있어 브라우저에서
  **동일한 파이썬 코드**를 실행함. 서버 버전과 결과가 100% 일치.
- DeepPrime 랭킹은 이 파일에는 포함되지 않음(genet은 파이썬 서버 전용). 체크박스를
  누르면 "서버 버전에서만 가능"이라는 안내가 뜨고, standalone 후보 목록은 그대로 나온다.

> 미리보기(iframe) 안에서 열면 `wasm instantiation failed` 같은 경고가 날 수 있는데,
> 이는 샌드박스 iframe이 WASM을 막아서 나는 것으로 **실제 브라우저 탭**에서 파일을
> 직접 열면 정상 동작한다.

### 2) 서버로 실행하기

> **요구사항: Python 3.10.** DeepPrime의 genet이 tensorflow<2.10을 요구해
> 3.11 이상에서는 설치되지 않는다. 서버만 쓸 때도 3.10으로 통일하는 것을 권장.

```bash
pip install -r requirements.txt      # 코어는 Flask만 필요
python run.py                        # http://127.0.0.1:8502 자동 오픈
```

`run.py`가 Flask 서버를 띄우고 브라우저를 연다. 수동으로 실행하려면 리포 루트에서:

```bash
PYTHONPATH=. python silent_mutation/webtool/server.py   # :8502
```

DeepPrime 효율 랭킹을 함께 사용하려면 (Python 3.10 환경에서) genet을 추가로 설치한다:

```bash
pip install -r requirements-deepprime.txt
```

없어도 웹툴은 정상 동작하며, DeepPrime만 "unavailable"로 표시된다.

---

## 디렉토리 구조

```
silent_mutation_standalone.html   설치 없이 여는 단일 파일 (1번 사용법)
run.py                            서버 실행 (2번 사용법)
build_standalone.py               단일 HTML 생성 스크립트 (코어 수정 시 재실행)
requirements.txt                  코어/서버 (Flask)
requirements-deepprime.txt        선택: DeepPrime (pandas + genet)
data/reference/codon_table.csv    코돈 테이블 (64개)
silent_mutation/
  core/     types, codon_utils, pam_finder, silent_finder, pegrna_builder, verify
  io/       sequence_loader (WT/Edited → Variant), genome_loader (reverse_complement),
            deepprime_runner (genet 효율 예측 — 서버 전용)
  webtool/  server.py (Flask), api.py (Flask-free 코어), index.html (서버 UI)
```

- **`api.py`**: `/api/analyze`·`/api/verify`의 로직을 Flask 없이 순수 함수
  (`run_analyze` / `run_verify`)로 뽑은 것. standalone HTML이 이걸 호출하며,
  server.py의 라우트와 **출력이 byte 단위로 동일함을 검증**했다.
- 코어(standalone 경로)는 pandas를 쓰지 않는다. codon 테이블을 stdlib `csv`로 읽는다.

---

## DeepPrime 연결

이 프로젝트의 범위는 **silent mutation 분석까지**이며, DeepPrime은 그 뒤에 붙는 별도
단계로 사용하였음. 추후 **deepPrime-coedit 툴**을 만들게 되면, 그때
변경해야 할 부분은 `silent_mutation/io/deepprime_runner.py`의 `run_deepprime_silent(...)` 이다.

  - 입력: `(wt, ed, ...pe_system, cell_type, rtt_max, top_n, exon_* )`
  - 출력: 효율 내림차순으로 정렬된 pegRNA dict 리스트. 각 dict는 그 pegRNA의 silent
    후보들을 `"outputs"`(list[PegRNAOutput])로 들고 있다.
  - 이 **입출력 형식만 유지**하면 내부 스코어링 엔진(genet → coedit 툴)만 갈아끼우면
    되며, `server.py`와 `index.html`은 변경할 필요 없음.

  동작 참고:
  - genet 없이 DeepPrime 랭킹을 요청하면,
    "DeepPrime 사용 불가(genet 설치 필요)"로 응답한다.
  - standalone HTML에는 DeepPrime이 포함되지 않는다(genet 미포함). 
    DeepPrime은 서버 빌드에서만 제공한다.

standalone 단일 파일은 `silent_mutation/` 코어의 **스냅샷**이다. 코어를 수정한 뒤에는
`python build_standalone.py`로 단일 파일을 다시 생성하면 된다 (코어 + `data/`를 gzip+base64로
HTML에 임베드하고, `fetch('/api/...')` 호출을 브라우저 내 `api.run_analyze`/`run_verify`
호출로 바꾼다).

---

## 설계 규약

- **Window 표준**: SynDesign/DeepPrime 표준 — 60bp flank, `VAR_IDX=60`, 1bp 변이는
  121bp 창.
- **좌표 규약**: `seq_wt`는 CDS 가닥 5'→3'; `pbs_seq`/`rtt_seq`는 PAM 가닥 5'→3'로 표기하도록 함;
  pegRNA 3' extension = `reverse_complement(pbs_seq + edited_rtt)`.
  
## 향후 확장 (미구현)

- **ClinVar ID 직접 입력**: 현재는 WT/Edited 서열 직접 입력만 지원한다. ClinVar ID로
  서열을 자동 조회하는 입력 모드는 아직 구현되지 않았다.
- **유전체 기반 배치 처리**: 여러 변이를 한 번에 처리하고 변이의 splice 경계
  (near_exon_boundary 등)를 자동 판정하는 파이프라인은 미구현이다. 이 기능은
  transcript/유전체 정보를 필요로 한다.

## 검증된 케이스

- EYS c.2528G>A (G843E, 121 candidates), EYS c.4957dupA (Ser1653 frameshift,
  37 candidates) — 유전체 EYS에서 추출해 end-to-end 검증.
- 서버 `/api/analyze`·`/api/verify` ↔ `api.run_analyze`/`run_verify` 출력 동일성 대조
  (DEMO 129 · BRCA2 27 · verify 포함) 통과.
