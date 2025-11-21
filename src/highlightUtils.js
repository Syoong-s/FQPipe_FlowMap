// src/highlightUtils.js

import Prism from "prismjs";
import "prismjs/components/prism-fortran";
import "prismjs/components/prism-python";

import functionCodeMap from "./subroutineCodeMap";

// 所有可跳转的函数名（来自函数实现表）
const functionNames = Object.keys(functionCodeMap);

// 生成带语法高亮的 HTML
export const getHighlightedHtml = (rawCode, language = "fortran") => {
  if (!rawCode) return "";
  const grammar = Prism.languages[language] || Prism.languages.plaintext;
  return Prism.highlight(rawCode, grammar, language);
};

// 在高亮后的 HTML 中，把函数名替换为可点击 span
export const addFunctionLinks = (html) => {
  if (!html) return "";
  let result = html;

  functionNames.forEach((name) => {
    const pattern = new RegExp(`\\b${name}\\b`, "g");
    result = result.replace(
      pattern,
      `<span class="func-link" data-func="${name}">${name}</span>`
    );
  });

  return result;
};
