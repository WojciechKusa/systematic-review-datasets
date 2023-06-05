import openai
import tiktoken

API_KEY_FILE = "../config/openai_key.txt"
with open(API_KEY_FILE, "r", encoding="utf-8") as f:
    api_key = f.read().strip()
openai.api_key = api_key



models = ["gpt-4-0314", "gpt-3.5-turbo-0301"]

model = models[1]
response = openai.ChatCompletion.create(
    model=model,
    messages=[
        {"role": "user", "content": "Does the following scientific paper fulfill all eligibility criteria and should it be included in the systematic review? Answer 'Included' or 'Excluded'."},
    ],
)

response

print(response)# x
