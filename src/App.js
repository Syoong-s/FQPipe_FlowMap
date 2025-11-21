import { useState, useEffect, useRef } from "react";
import ReactFlow, { Background, useReactFlow, ReactFlowProvider } from "reactflow";
import "reactflow/dist/style.css";
// ==========引入配色=========================
import Prism from "prismjs";
import "prismjs/components/prism-fortran";
import "prismjs/components/prism-python";
import "./styles/prism_vscode.css";
// ==========引入初始代码=======================
import initialCode from "./initialCode";
import codeMapMain from "./codeMapMain";
import codeMapSub from "./codeMapSub";

// 公共样式：让每个节点大一点、字体大一点
const nodeStyle = {
  width: 260,
  height: 70,
  fontSize: 20,
  fontWeight: "bold",
  display: "flex",
  justifyContent: "center",
  alignItems: "center",
   border: "1px solid #555",
  borderRadius: 8,
  background: "#2d2d2d",    // 节点暗色
  color: "#eaeaea",         // 白色字体
   borderRadius: 6
};

const nodes = [
  { id: 'chip_pre', position: { x: -400, y: -450 }, data: { label: "pre_process" } ,style: nodeStyle },
  { id: 'astro', position: { x: -400, y: -300 }, data: { label: "proc_astrometry" } ,style: nodeStyle },
  { id: 'source', position: { x: -400, y: -150 }, data: { label: "proc_source" } ,style: nodeStyle },
  { id: 'fourier', position: { x: -400, y: 0 }, data: { label: "proc_Fourier_T" } ,style: nodeStyle },
  { id: 'psf', position: { x: -400, y: 150 }, data: { label: "PSF_reconstruction" } ,style: nodeStyle },
  { id: 'shear', position: { x: -400, y: 300 }, data: { label: "gen_shear_cat" } ,style: nodeStyle },
  { id: 'combine', position: { x: -400, y: 450 }, data: { label: "combine_expo_catalog" } ,style: nodeStyle }
];
const mainTitleMap = {
  chip_pre: "pre_process",
  astro: "proc_astrometry",
  source: "proc_source",
  fourier: "proc_Fourier_T",
  psf: "PSF_reconstruction",
  shear: "gen_shear_cat",
  combine: "combine_expo_catalog"
};
const edges = [
  { id: 'e1', source: 'chip_pre', target: 'astro' },
  { id: 'e2', source: 'astro', target: 'source' },
  { id: 'e3', source: 'source', target: 'fourier' },
  { id: 'e4', source: 'fourier', target: 'psf' },
  { id: 'e5', source: 'psf', target: 'shear' },
  { id: 'e6', source: 'shear', target: 'combine' }
];

const subNodeStyle = {
  width: 260,
  height: 55,
  fontSize: 18,
  fontWeight: "bold",
  display: "flex",
  justifyContent: "center",
  alignItems: "center",
  background: "#2d2d2d",
  color: "#eaeaea",
  border: "1px solid #666",
  borderRadius: 6
};

const subFlowMap = {

  chip_pre: {
    nodes: [
      { id: "Set_Background", position: { x: 0, y: -150 }, data:{ label: "Set_Background" }, style: subNodeStyle },
      { id: "Star_match", position: { x: 0, y: -60 }, data:{ label: "Star_match" }, style: subNodeStyle },
      { id: "Locate_defects", position: { x: 0, y: 30 }, data:{ label: "Locate_defects" }, style: subNodeStyle },
      { id: "Merge_defects", position: { x: 0, y: 120 }, data:{ label: "Merge_defects" }, style: subNodeStyle },
      { id: "Gen_norm", position: { x: 0, y: 210 }, data:{ label: "Gen_norm.fits" }, style: subNodeStyle },
    ],
    edges: [
      { id: "p1", source: "Set_Background", target: "Star_match" },
      { id: "p2", source: "Star_match", target: "Locate_defects" },
      { id: "p3", source: "Locate_defects", target: "Merge_defects" },
      { id: "p4", source: "Merge_defects", target: "Gen_norm" }
    ]
  },

  // ================= proc_astrometry =================
  astro: {
    nodes: [
      { id: "get_astrometry", position: { x: 0, y: -80 }, data:{ label: "get_astrometry" }, style: subNodeStyle },
      { id: "update_norm", position: { x: 0, y: 20 }, data:{ label: "update_norm" }, style: subNodeStyle },
    ],
    edges: [
      { id: "a1", source: "get_astrometry", target: "update_norm" }
    ]
  },

  // ================= proc_source =================
  source: {
    nodes: [
      { id: "get_weight", position: { x: 0, y: -150 }, data:{ label: "get_weight" }, style: subNodeStyle },
      { id: "get_expo_catalog", position: { x: 0, y: -60 }, data:{ label: "get_expo_catalog" }, style: subNodeStyle },
      { id: "gen_source_ext_catalog", position: { x: 0, y: 30 }, data:{ label: "gen_source_ext_catalog" }, style: subNodeStyle },
      { id: "gen_star_candidate_direct", position: { x: 0, y: 120 }, data:{ label: "gen_star_candidate_direct" }, style: subNodeStyle },
    ],
    edges: [
      { id: "s1", source: "get_weight", target: "get_expo_catalog" },
      { id: "s2", source: "get_expo_catalog", target: "gen_source_ext_catalog" },
      { id: "s3", source: "gen_source_ext_catalog", target: "gen_star_candidate_direct" }
    ]
  },

  // ================= proc_Fourier_T =================
  fourier: {
    nodes: [
      { id: "get_gal_power", position: { x: 0, y: -60 }, data:{ label: "get_gal_power" }, style: subNodeStyle },
      { id: "get_star_power", position: { x: 0, y: 40 }, data:{ label: "get_star_power" }, style: subNodeStyle },
    ],
    edges: [
      { id: "f1", source: "get_gal_power", target: "get_star_power" }
    ]
  },

  // ================= PSF_reconstruction =================
  psf: {
    nodes: [
      { id: "read_in_candidates", position: { x: 0, y: -180 }, data:{ label: "read_in_candidates" }, style: subNodeStyle },
      { id: "star_selection", position: { x: 0, y: -90 }, data:{ label: "star_selection" }, style: subNodeStyle },
      { id: "plot_star_expo", position: { x: 0, y: 0 }, data:{ label: "plot_star_expo" }, style: subNodeStyle },
      { id: "plot_stars", position: { x: 0, y: 90 }, data:{ label: "plot_stars" }, style: subNodeStyle },
      { id: "make_PSF_local_fit", position: { x: 0, y: 180 }, data:{ label: "make_PSF_local_fit" }, style: subNodeStyle },
    ],
    edges: [
      { id: "ps1", source: "read_in_candidates", target: "star_selection" },
      { id: "ps2", source: "star_selection", target: "plot_star_expo" },
      { id: "ps3", source: "plot_star_expo", target: "plot_stars" },
      { id: "ps4", source: "plot_stars", target: "make_PSF_local_fit" }
    ]
  },

  // ================= gen_shear_cat =================
  shear: {
    nodes: [
      { id: "Load_PSF", position: { x: 0, y: -150 }, data:{ label: "Load_PSF" }, style: subNodeStyle },
      { id: "read_astrometry_para", position: { x: 0, y: -60 }, data:{ label: "read_astrometry_para" }, style: subNodeStyle },
      { id: "read_source_info", position: { x: 0, y: 30 }, data:{ label: "read_source_info" }, style: subNodeStyle },
      { id: "field_distortion_PU", position: { x: 0, y: 120 }, data:{ label: "field_distortion_PU" }, style: subNodeStyle },
      { id: "get_shear", position: { x: 0, y: 210 }, data:{ label: "get_shear" }, style: subNodeStyle },
    ],
    edges: [
      { id: "sh1", source: "Load_PSF", target: "read_astrometry_para" },
      { id: "sh2", source: "read_astrometry_para", target: "read_source_info" },
      { id: "sh3", source: "read_source_info", target: "field_distortion_PU" },
      { id: "sh4", source: "field_distortion_PU", target: "get_shear" }
    ]
  },

  // ================= combine_expo_catalog =================
  combine: {
    nodes: [
      { id: "read_sheardat", position: { x: 0, y: -60 }, data: { label: "read_sheardat" }, style: subNodeStyle },
      { id: "Filtering_Combining", position: { x: 0, y: 40 }, data:{ label: "Filtering_Combining" }, style: subNodeStyle },
    ],
    edges: [
      { id: "c1", source: "read_sheardat", target: "Filtering_Combining" }
    ]
  }
};

// ================== 组件主体 ==================

function MainApp() {
  // 状态管理
  const [codeHtml, setCodeHtml] = useState(""); // 这里存的是处理好的 HTML
  const [activeMain, setActiveMain] = useState(null);
  const codeRef = useRef(null);

  // ★★★ 辅助函数：生成带语法高亮的 HTML ★★★
  const getHighlightedHtml = (rawCode, language = "fortran") => {
    if (!rawCode) return "";
    // 确保 Prism 语言包已加载，如果没有加载 fortran，回退到 plain text
    const grammar = Prism.languages[language] || Prism.languages.plaintext;
    return Prism.highlight(rawCode, grammar, language);
  };

  // 初始化：组件挂载时高亮 initialCode
  useEffect(() => {
    setCodeHtml(getHighlightedHtml(initialCode));
  }, []);

  // ★★★ 修复后的子模块高亮逻辑 ★★★
  function highlightSubmodule(subId) {
    console.log("Clicked submodule:", subId);
    
    // 1. 获取范围
    const range = codeMapSub[subId];
    if (!range) {
      console.warn(`未找到子模块 ${subId} 的代码映射范围`);
      return;
    }
    
    // 2. 获取当前主模块的完整源代码
    const rawCode = codeMapMain[activeMain];
    if (!rawCode) return;

    // 3. 先进行 Prism 语法高亮，得到带 <span> 颜色的 HTML 字符串
    const syntaxHtml = getHighlightedHtml(rawCode, "fortran");

    // 4. 按行拆分
    const lines = syntaxHtml.split("\n");

    // 5. 包裹高亮区域
    for (let i = range.start - 1; i <= range.end - 1; i++) {
      if (lines[i] !== undefined) {
        // 使用 display: inline-block 确保背景色占满整行
        lines[i] = `<mark class="highlight-block" style="background: #5c4e22; color: inherit; display: inline-block; width: 100%; border-left: 3px solid #f1c40f; padding-left: 4px;">${lines[i]}</mark>`;
      }
    }

    // 6. 合并回字符串并更新状态
    const finalHtml = lines.join("\n");
    setCodeHtml(finalHtml);

    // 7. 滚动逻辑
    setTimeout(() => {
      const mark = document.querySelector("mark.highlight-block");
      if (mark) {
        mark.scrollIntoView({ behavior: "smooth", block: "center" });
      }
    }, 50);
  }

  return (
    <div style={{ display: "flex", height: "100vh", overflow: "hidden", background: "#1e1e1e", color: "#eaeaea" }}>

      {/* 左侧流程图区域 */}
      <div style={{
        flex: 25, backgroundColor: "#1e1e1e", borderLeft: "2px solid #444",
        borderRight: "2px solid #444", position: "relative",
        display: "flex", flexDirection: "column",
         minWidth: 0
      }}>
        <div style={{
          padding: "15px 30px", fontSize: "28px", fontWeight: "bold", color: "#eaeaea",
          borderBottom: "1px solid #444", letterSpacing: "2px",
        }}>
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
                // ★ 点击主节点，生成完整代码的高亮HTML
                const rawCode = codeMapMain[node.id] || "";
                setCodeHtml(getHighlightedHtml(rawCode, "fortran"));
              }}
              proOptions={{ hideAttribution: true }}
            >
              <Background color="#aaa" gap={16} />
            </ReactFlow>
          </ReactFlowProvider>
        </div>
      </div>

      {/* 中间子模块流程图 */}
      <div style={{
        flex: 20, backgroundColor: "#1e1e1e", borderLeft: "2px solid #444",
        borderRight: "2px solid #444", overflow: "hidden", position: "relative",
        display: "flex", flexDirection: "column",
         minWidth: 0
      }}>
        <div style={{
          padding: "15px 30px", fontSize: "28px", fontWeight: "bold", color: "#eaeaea",
          borderBottom: "1px solid #444", letterSpacing: "2px",
        }}>
          Subroutine
        </div>

        <div style={{ flex: 1 }}>
          <ReactFlowProvider>
            <ReactFlow
              id="subFlow"
              nodes={activeMain && subFlowMap[activeMain] ? subFlowMap[activeMain].nodes : []}
              edges={activeMain && subFlowMap[activeMain] ? subFlowMap[activeMain].edges : []}
              onNodeClick={(evt, node) => {
                evt.stopPropagation();
                // ★ 点击子节点，触发局部高亮
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

      {/* 右侧代码区 */}
      <div style={{
        flex: 50, display: "flex", flexDirection: "column",
        background: "#1e1e1e", borderLeft: "2px solid #8d8484ff",
        minWidth: 0 // ★ 防止内容撑开父容器
      }}>
        <div style={{
          padding: "15px 30px", fontSize: "32px", fontWeight: "bold", color: "#eaeaea",
          borderBottom: "1px solid #444", letterSpacing: "2px",
        }}>
          {activeMain ? mainTitleMap[activeMain] : "MainPipeline"}
        </div>

        <pre
          ref={codeRef}
          style={{
            flex: 1, margin: 0, padding: 40, background: "#1e1e1e",
            overflow: "auto", whiteSpace: "pre", // 不换行
            color: "#f7f2f2ff", boxSizing: "border-box",
            fontFamily: 'Consolas, Monaco, "Andale Mono", "Ubuntu Mono", monospace',
            fontSize: '14px',
            lineHeight: '1.5'              
          }}
        >
          <code
            className="language-fortran"
            dangerouslySetInnerHTML={{
              __html: codeHtml || `<span style="color: #888;">点击左侧主模块或子模块以显示代码...</span>`,
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