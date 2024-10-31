Shader "Unlit/RaymarchMetaball"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _Color ("Color", Color) = (1, 1, 1, 1)
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" "Queue"="Transparent" }
        Blend SrcAlpha OneMinusSrcAlpha
        ZTest Always
        ZWrite OFF
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            #define MAX_STEPS 100
            #define MAX_DIST 100
            #define SURF_DIST 0.001

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
                float3 ro : TEXCOORD1;
                float3 hitPos : TEXCOORD2;
                float3 light : TEXCOORD3;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;
            sampler2D _CameraDepthTexture;
            float4 _Color;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                o.ro = mul(unity_WorldToObject, float4(_WorldSpaceCameraPos, 1));
                o.hitPos = v.vertex;
                o.light = mul(unity_WorldToObject, _WorldSpaceLightPos0);
                return o;
            }

            float Metaball(float3 d, float3 c, float r) {
                return r / (length(d-c));   
            }

            // exponential
            float smin( float a, float b, float c, float k )
            {
                k *= 1.0;
                float r = exp2(-a/k) + exp2(-b/k) + exp2(-c/k);
                return -k*log2(r);
            }

            float GetDist(float3 p) {
                float3 o1 = float3(0.4, 0, 0);
                float3 o2 = float3(-0.4, 0, 0);
                float3 o3 = float3(0, 0.3, 0);
                
                float d1 = length(p - o1) - 0.2; // sphere 1
                float d2 = length(p - o2) - 0.2; // sphere 2
                float d3 = length(p - o3) - 0.2;
                float d = smin(d1, d2, d3, 0.1);

                //d = -Metaball(p, float3(0, 0, 0), 0.5) + 1;
                //d += -Metaball(p, float3(0.2, 0, 0), 0.5) + 1;

                return d;
            }

            float Raymarch(float3 ro, float3 rd) {
                float dO = 0;
                float dS;
                for (int i = 0; i < MAX_STEPS; i++) {
                    float3 p = ro + dO * rd;
                    dS = GetDist(p);
                    dO += dS;

                    if (dS < SURF_DIST || dO > MAX_DIST) break;
                }
                return dO;
            }

            float3 GetNormal(float3 p) {
                float2 e = float2(1e-2, 0);
                float3 n = GetDist(p) - float3(
                    GetDist(p - e.xyy),
                    GetDist(p - e.yxy),
                    GetDist(p - e.yyx)
                    );
                return normalize(n);
            }

            fixed4 frag (v2f i) : SV_Target
            {
                float2 uv = i.uv - 0.5;
                float3 ro = i.ro;
                float3 rd = normalize(i.hitPos - ro);

                float d = Raymarch(ro, rd);
                fixed4 tex = tex2D(_MainTex, i.uv);
                fixed4 col = 0;
                float m = dot(uv, uv);

                if (d >= MAX_DIST) {
                    discard;
                } else {
                    float3 p = ro + rd * d;
                    float3 n = GetNormal(p);
                    float cosNR = dot(n, normalize(ro));
                    float sinNR = pow(sqrt(1 - cosNR * cosNR), 3);
                    col.w = 1;
                    col.rgba = _Color.rgba * sinNR;
                    
                    // depth test
                    float4 world = mul(unity_ObjectToWorld, float4(p, 1));
                    float4 view = mul(UNITY_MATRIX_V, world);
                    float4 proj = mul(UNITY_MATRIX_P, view);
                    proj /= proj.w;

                    float sceneDepth = tex2D(_CameraDepthTexture, float2(0.5 + proj.x / 2, 0.5 - proj.y / 2));
                    if (proj.z < sceneDepth)
                        discard;

                    // specular
                    float3 h = -normalize(rd - i.light);
                    float sp = pow(saturate(dot(h, n)), 20);
                    col += sp * 0.7;
                }
                return col;
            }

            
            // https://www.youtube.com/watch?v=S8AWd66hoCo&ab_channel=TheArtofCode
            ENDCG
        }
    }
}
