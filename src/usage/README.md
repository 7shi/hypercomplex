# usage

LLM呼び出しのトークン使用量（[llm7shi](https://github.com/7shi/llm7shi)の`Usage`）をリポジトリ直下の[usage.txt](../../usage.txt)に記録し、日付ごとに集計するツールです。

```bash
uv run usage [-a]
```

| 引数 | 説明 | デフォルト |
|---|---|---|
| `-a`, `--all` | 日付ごとの合計をすべて表示し、最後に全期間の総合計を表示 | 今日の分のみ表示 |

## usage.txtの形式

`# YYYY/MM/DD`の見出しと、その日に実行した各コマンドの`Usage`の`repr`を並べたテキストファイルです。

```
# 2026/09/16
Usage({'input_tokens': 17221, 'output_tokens': 2225, 'reasoning_tokens': 169, 'cached_tokens': 0, 'total_tokens': 19446})
Usage({'input_tokens': 27100, 'output_tokens': 2688, 'reasoning_tokens': 447, 'cached_tokens': 0, 'total_tokens': 29788})

# 2026/09/17
Usage({'input_tokens': 24402, 'output_tokens': 1830, 'reasoning_tokens': 196, 'cached_tokens': 0, 'total_tokens': 26232})
```

## API

他のツールから使う関数です。

| 関数 | 説明 |
|---|---|
| `parse_usage_file(path=USAGE_PATH)` | ファイルをパースし、日付をkey、その日の合計`Usage`をvalueとするdictを返す（日付は出現順） |
| `append_usage(usage, path=USAGE_PATH, date=None)` | `Usage`を指定日（省略時は今日）のセクションに追記する。見出しがなければ末尾に作成する |

[review](../review/README.md)は実行の最後に`append_usage`でその実行の合計を追記します。
