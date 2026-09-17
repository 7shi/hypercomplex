# usage

LLM呼び出しのトークン使用量（[llm7shi](https://github.com/7shi/llm7shi)の`Usage`）をリポジトリ直下の[usage.jsonl](../../usage.jsonl)に記録し、日付ごとに集計するツールです。

```bash
uv run usage [-a]
```

| 引数 | 説明 | デフォルト |
|---|---|---|
| `-a`, `--all` | 日付ごとの合計をすべて表示し、最後に全期間の総合計を表示 | 今日の分のみ表示 |

## usage.jsonlの形式

1行1レコードのJSON Lines形式です。各レコードはタイムゾーン付きの生成日時（`timestamp`）、モデル名（`model`）、`Usage`の`to_dict()`（`raw`を除いた正規化済みフィールド）を持ちます。

```
{"timestamp": "2026-09-16T00:00:00+00:00", "model": "gpt-6-astra", "input_tokens": 17221, "output_tokens": 2225, "reasoning_tokens": 169, "cached_tokens": 0, "total_tokens": 19446}
{"timestamp": "2026-09-17T09:12:34+09:00", "model": "gpt-6-astra", "input_tokens": 24402, "output_tokens": 1830, "reasoning_tokens": 196, "cached_tokens": 0, "total_tokens": 26232}
```

記録はローカルのタイムゾーン付きで行うが、日付ごとの集計は各レコードの`timestamp`をUTCに変換してから日付部分を使う。既存データ移行時など生成時刻が不明なレコードはUTC（`+00:00`）の00:00:00として記録している。

集計は日付×モデルの二軸で行い、`--all`指定時は末尾にモデルごとの全期間合計も表示する。

```
2026/09/16
  gpt-6-astra {'input_tokens': ...}
2026/09/17
  gpt-6-astra {'input_tokens': ...}
==========
gpt-6-astra {'input_tokens': ...}
```

## API

他のツールから使う関数です。

| 関数 | 説明 |
|---|---|
| `parse_usage_file(path=USAGE_PATH)` | ファイルをパースし、日付をkey、モデル名をkeyとする合計`Usage`のdictをvalueとするdictを返す（日付・モデル名は出現順） |
| `append_usage(usage, model, path=USAGE_PATH, timestamp=None)` | `Usage`をモデル名・タイムゾーン付きの生成日時（省略時は現在時刻）とともに1行追記する |

[review](../review/README.md)は実行の最後に`append_usage`でその実行の合計を追記します。
