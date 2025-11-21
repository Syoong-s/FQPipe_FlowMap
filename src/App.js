// src/App.js

import { useState, useEffect, useRef } from "react";
import ReactFlow, { Background, ReactFlowProvider } from "reactflow";
import "reactflow/dist/style.css";

import "./styles/prism_vscode.css";

import initialCode from "./initialCode";
import codeMapMain from "./codeMapMain";
import codeMapSub from "./codeMapSub";
import functionCodeMap from "./subroutineCodeMap";

import { nodes, edges, subFlowMap, mainTitleMap } from "./flowConfig";
import { getHighlightedHtml, addFunctionLinks } from "./highlightUtils";

// ================== 组件主体 ==================

function MainApp() {
  const [codeHtml, setCodeHtml] = useState("");
  const [activeMain, setActiveMain] = useState(null);
  const codeRef = useRef(null);

  // 右上角标题 & 代码视图栈（用于返回）
  const [viewStack, setViewStack] = useState([]);
  const [currentTitle, setCurrentTitle] = useState("MainPipeline");

  // 初始化：组件挂载时高亮 initialCode
  useEffect(() => {
    const html = getHighlightedHtml(initialCode, "fortran");
    setCodeHtml(addFunctionLinks(html));
    setCurrentTitle("MainPipeline");
  }, []);

  // 点击子流程节点 → 高亮主代码里的对应行
  function highlightSubmodule(subId) {
    console.log("Clicked submodule:", subId);

    const range = codeMapSub[subId];
    if (!range) {
      console.warn(`未找到子模块 ${subId} 的代码映射范围`);
      return;
    }

    const rawCode = codeMapMain[activeMain];
    if (!rawCode) return;

    const syntaxHtml = getHighlightedHtml(rawCode, "fortran");
    const lines = syntaxHtml.split("\n");

    for (let i = range.start - 1; i <= range.end - 1; i++) {
      if (lines[i] !== undefined) {
        lines[i] = `<mark class="highlight-block" style="background: #5c4e22; color: inherit; display: inline-block; width: 100%; border-left: 3px solid #f1c40f; padding-left: 4px;">${lines[i]}</mark>`;
      }
    }

    const finalHtml = lines.join("\n");
    const finalWithLinks = addFunctionLinks(finalHtml);
    setCodeHtml(finalWithLinks);

    setTimeout(() => {
      const mark = document.querySelector("mark.highlight-block");
      if (mark) {
        mark.scrollIntoView({ behavior: "smooth", block: "center" });
      }
    }, 50);
  }

  // Ctrl+点击函数名 → 跳转到函数实现
  function handleFunctionJump(funcName) {
    const implCode = functionCodeMap[funcName];
    if (!implCode) {
      console.warn("未找到函数实现：", funcName);
      return;
    }

    const currentScrollTop = codeRef.current ? codeRef.current.scrollTop : 0;

    setViewStack((prev) => [
      ...prev,
      {
        html: codeHtml,
        title: currentTitle,
        scrollTop: currentScrollTop,
      },
    ]);

    setCurrentTitle(funcName);

    const syntaxHtml = getHighlightedHtml(implCode, "fortran");
    const withLinks = addFunctionLinks(syntaxHtml);
    setCodeHtml(withLinks);

    setTimeout(() => {
      if (codeRef.current) {
        codeRef.current.scrollTop = 0;
      }
    }, 0);
  }

  // 为代码区绑定 Ctrl+点击事件
  useEffect(() => {
    const container = codeRef.current;
    if (!container) return;

    const handleClick = (e) => {
      const isJumpKey = e.ctrlKey || e.metaKey;
      if (!isJumpKey) return;

      let el = e.target;
      while (el && el !== container) {
        if (el.getAttribute && el.getAttribute("data-func")) {
          const funcName = el.getAttribute("data-func");
          e.preventDefault();
          handleFunctionJump(funcName);
          break;
        }
        el = el.parentElement;
      }
    };

    container.addEventListener("click", handleClick);
    return () => container.removeEventListener("click", handleClick);
  }, [codeHtml, currentTitle]);

  return (
    <div
      style={{
        display: "flex",
        height: "100vh",
        overflow: "hidden",
        background: "#1e1e1e",
        color: "#eaeaea",
      }}
    >
      {/* 左侧：主流程图 */}
      <div
        style={{
          flex: 25,
          backgroundColor: "#1e1e1e",
          borderLeft: "2px solid #444",
          borderRight: "2px solid #444",
          position: "relative",
          display: "flex",
          flexDirection: "column",
          minWidth: 0,
        }}
      >
        <div
          style={{
            padding: "15px 30px",
            fontSize: "28px",
            fontWeight: "bold",
            color: "#eaeaea",
            borderBottom: "1px solid #444",
            letterSpacing: "2px",
          }}
        >
          MainPipe
        </div>

        <div style={{ flex: 1 }}>
          <ReactFlowProvider>
            <ReactFlow
              id="mainFlow"
              nodes={nodes}
              edges={edges}
              nodesDraggable={false}
              nodesConnectable={false}
              zoomOnScroll={false}
              zoomOnPinch={false}
              zoomOnDoubleClick={false}
              panOnDrag={false}
              panOnScroll={false}
              panOnZoom={false}
              fitView
              onNodeClick={(evt, node) => {
                evt.stopPropagation();
                setActiveMain(node.id);

                setViewStack([]);
                setCurrentTitle(mainTitleMap[node.id] || "MainPipeline");

                const rawCode = codeMapMain[node.id] || "";
                const syntaxHtml = getHighlightedHtml(rawCode, "fortran");
                setCodeHtml(addFunctionLinks(syntaxHtml));
              }}
              proOptions={{ hideAttribution: true }}
            >
              <Background color="#aaa" gap={16} />
            </ReactFlow>
          </ReactFlowProvider>
        </div>
      </div>

      {/* 中间：子流程图 */}
      <div
        style={{
          flex: 20,
          backgroundColor: "#1e1e1e",
          borderLeft: "2px solid #444",
          borderRight: "2px solid #444",
          overflow: "hidden",
          position: "relative",
          display: "flex",
          flexDirection: "column",
          minWidth: 0,
        }}
      >
        <div
          style={{
            padding: "15px 30px",
            fontSize: "28px",
            fontWeight: "bold",
            color: "#eaeaea",
            borderBottom: "1px solid #444",
            letterSpacing: "2px",
          }}
        >
          Subprocess
        </div>

        <div style={{ flex: 1 }}>
          <ReactFlowProvider>
            <ReactFlow
              id="subFlow"
              nodes={
                activeMain && subFlowMap[activeMain]
                  ? subFlowMap[activeMain].nodes
                  : []
              }
              edges={
                activeMain && subFlowMap[activeMain]
                  ? subFlowMap[activeMain].edges
                  : []
              }
              onNodeClick={(evt, node) => {
                evt.stopPropagation();
                highlightSubmodule(node.id);
              }}
              nodesDraggable={false}
              nodesConnectable={false}
              zoomOnScroll={false}
              zoomOnPinch={false}
              zoomOnDoubleClick={false}
              panOnDrag={false}
              panOnScroll={false}
              panOnZoom={false}
              fitView
              proOptions={{ hideAttribution: true }}
            >
              <Background color="#555" gap={16} />
            </ReactFlow>
          </ReactFlowProvider>
        </div>
      </div>

      {/* 右侧：代码区 */}
      <div
        style={{
          flex: 50,
          display: "flex",
          flexDirection: "column",
          background: "#1e1e1e",
          borderLeft: "2px solid #8d8484ff",
          minWidth: 0,
        }}
      >
        <div
          style={{
            padding: "15px 30px",
            fontSize: "32px",
            fontWeight: "bold",
            color: "#eaeaea",
            borderBottom: "1px solid #444",
            letterSpacing: "2px",
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
          }}
        >
          {/* 标题 */}
          <span>{currentTitle}</span>

          {/* 右侧按钮区：返回 + 初始代码 */}
          <div style={{ display: "flex", gap: "8px" }}>
            {viewStack.length > 0 && (
              <button
                onClick={() => {
                  setViewStack((prev) => {
                    if (prev.length === 0) return prev;
                    const last = prev[prev.length - 1];

                    setCurrentTitle(last.title);
                    setCodeHtml(last.html);

                    setTimeout(() => {
                      if (codeRef.current) {
                        codeRef.current.scrollTop = last.scrollTop;
                      }
                    }, 0);

                    return prev.slice(0, -1);
                  });
                }}
                style={{
                  padding: "4px 12px",
                  fontSize: "14px",
                  background: "#333",
                  color: "#eee",
                  borderRadius: 4,
                  border: "1px solid #666",
                  cursor: "pointer",
                }}
              >
                返回
              </button>
            )}

            <button
              onClick={() => {
                setViewStack([]);
                setCurrentTitle("MainPipeline");

                const html = getHighlightedHtml(initialCode, "fortran");
                const withLinks = addFunctionLinks(html);
                setCodeHtml(withLinks);

                if (codeRef.current) {
                  codeRef.current.scrollTop = 0;
                }
              }}
              style={{
                padding: "4px 12px",
                fontSize: "14px",
                background: "#444",
                color: "#eee",
                borderRadius: 4,
                border: "1px solid #777",
                cursor: "pointer",
              }}
            >
              初始代码
            </button>
          </div>
        </div>

        <pre
          ref={codeRef}
          style={{
            flex: 1,
            margin: 0,
            padding: 40,
            background: "#1e1e1e",
            overflow: "auto",
            whiteSpace: "pre",
            color: "#f7f2f2ff",
            boxSizing: "border-box",
            fontFamily:
              'Consolas, Monaco, "Andale Mono", "Ubuntu Mono", monospace',
            fontSize: "14px",
            lineHeight: "1.5",
          }}
        >
          <code
            className="language-fortran"
            dangerouslySetInnerHTML={{
              __html:
                codeHtml ||
                `<span style="color: #888;">点击左侧主模块或子模块以显示代码...</span>`,
            }}
          />
        </pre>
      </div>
    </div>
  );
}

export default function App() {
  return <MainApp />;
}
